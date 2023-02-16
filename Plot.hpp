#pragma once

#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdarg>
#include <cassert>

#ifdef WIN32
#include <io.h>
#endif // #ifdef WIN32

#define GP_MAX_TMP_FILES    64


typedef struct _GNUPLOT_CTRL_ {
    /** Pipe to gnuplot process */
    FILE    * gnucmd ;

    /** Number of currently active plots */
    int       nplots ;
    /** Current plotting style */
    char      pstyle[32] ;

    /** Pointer to table of names of temporary files */
    char*      tmp_filename_tbl[GP_MAX_TMP_FILES] ;
    /** Number of temporary files */
    int       ntmp ;
} gnuplot_ctrl ;

template<typename T>
struct Plot
{
    gnuplot_ctrl *handle;

    Plot() {        
        int i;

    #ifndef WIN32
        if (getenv("DISPLAY") == NULL) {
            fprintf(stderr, "cannot find DISPLAY variable: is it set?\n") ;
        }
    #endif // #ifndef WIN32


        /*
        * Structure initialization:
        */
        handle = (gnuplot_ctrl*)malloc(sizeof(gnuplot_ctrl)) ;
        handle->nplots = 0 ;
        setstyle("points") ;
        handle->ntmp = 0 ;

        handle->gnucmd = popen("gnuplot 2> /dev/null", "w") ;
        if (handle->gnucmd == NULL) {
            fprintf(stderr, "error starting gnuplot, is gnuplot or gnuplot.exe in your path?\n") ;
            free(handle) ;
            exit(-1);
        }

        for (i=0;i<GP_MAX_TMP_FILES; i++)
        {
            handle->tmp_filename_tbl[i] = NULL;
        }        
    }


    ~Plot() {
       int     i ;
        if (pclose(handle->gnucmd) == -1) {
            fprintf(stderr, "problem closing communication to gnuplot\n") ;
            return ;
        }
        if (handle->ntmp) {
            for (i=0 ; i<handle->ntmp ; i++) {
                remove(handle->tmp_filename_tbl[i]) ;
                free(handle->tmp_filename_tbl[i]);
                handle->tmp_filename_tbl[i] = NULL;

            }
        }
        free(handle) ;
        return ;
    }

    void cmd(char const *  cmd, ...)
    {
        va_list ap ;
        va_start(ap, cmd);
        vfprintf(handle->gnucmd, cmd, ap);
        va_end(ap);

        fputs("\n", handle->gnucmd) ;
        fflush(handle->gnucmd) ;
        return ;
    }

    void setstyle(const char * plot_style)
    {
        if (strcmp(plot_style, "lines") &&
            strcmp(plot_style, "points") &&
            strcmp(plot_style, "linespoints") &&
            strcmp(plot_style, "impulses") &&
            strcmp(plot_style, "dots") &&
            strcmp(plot_style, "steps") &&
            strcmp(plot_style, "errorbars") &&
            strcmp(plot_style, "boxes") &&
            strcmp(plot_style, "boxerrorbars")) {
            fprintf(stderr, "warning: unknown requested style: using points\n") ;
            strcpy(handle->pstyle, "points") ;
        } else {
            strcpy(handle->pstyle, plot_style) ;
        }
        return ;
    }

    void set_xlabel(char * label)
    {
        cmd("set xlabel \"%s\"", label) ;
    }

    void set_ylabel(char * label)
    {
        cmd("set ylabel \"%s\"", label) ;
    }

    void resetplot()
    {
        int     i ;
        if (handle->ntmp) {
            for (i=0 ; i < handle->ntmp ; i++) {
                remove( handle->tmp_filename_tbl[i]) ;
                free(  handle->tmp_filename_tbl[i]);
                handle->tmp_filename_tbl[i] = NULL;

            }
        }
        handle->ntmp = 0 ;
        handle->nplots = 0 ;
        return ;
    }

    void plot_x(    
    T              *   d,
    int                 n,
    const char            *   title
    )
    {
        int     i ;
        FILE*   tmpfd ;
        char const * tmpfname;

        if (d==NULL || (n<1)) return ;

        /* Open temporary file for output   */
        tmpfname = tmpfile();
        tmpfd = fopen(tmpfname, "w");

        if (tmpfd == NULL) {
            fprintf(stderr,"cannot create temporary file: exiting plot") ;
            return ;
        }

        /* Write data to this file  */
        for (i=0 ; i<n ; i++) {
            fprintf(tmpfd, "%.18e\n", d[i]);
        }
        fclose(tmpfd) ;

        plot_atmpfile(tmpfname,title);
        return ;
    }

    void plot_xy(        
        T          *   x,
        T          *   y,
        int                 n,
        char            *   title
    )
    {
        int     i ;
        FILE*   tmpfd ;
        char const * tmpfname;

        if (x==NULL || y==NULL || (n<1)) return ;

        /* Open temporary file for output   */
        tmpfname = tmpfile();
        tmpfd = fopen(tmpfname, "w");

        if (tmpfd == NULL) {
            fprintf(stderr,"cannot create temporary file: exiting plot") ;
            return ;
        }

        /* Write data to this file  */
        for (i=0 ; i<n; i++) {
            fprintf(tmpfd, "%.18e %.18e\n", x[i], y[i]) ;
        }
        fclose(tmpfd) ;

        plot_atmpfile(tmpfname,title);
        return ;
    }

    void plot_once(
    char    *   title,
    char    *   style,
    char    *   label_x,
    char    *   label_y,
    T  *   x,
    T  *   y,
    int         n
    )
    {
    if (x==NULL || n<1) return ;    
    if (style!=NULL) {
        setstyle( style);
    } else {
        setstyle( "lines");
    }
    if (label_x!=NULL) {
        set_xlabel( label_x);
    } else {
        set_xlabel( "X");
    }
    if (label_y!=NULL) {
        set_ylabel( label_y);
    } else {
        set_ylabel( "Y");
    }
    if (y==NULL) {
        plot_x( x, n, title);
    } else {
        plot_xy( x, y, n, title);
    }
    printf("press ENTER to continue\n");
    while (getchar()!='\n') {}
    close(handle);
    return ;
    }

    void plot_slope(        
        T              a,
        T              b,
        char            *   title
    )
    {
        char const *    scmd    = (handle->nplots > 0) ? "replot" : "plot";
        title                  = (title == NULL)      ? "(none)" : title;

        cmd( "%s %.18e * x + %.18e title \"%s\" with %s",
                    scmd, a, b, title, handle->pstyle) ;

        handle->nplots++ ;
        return ;
    }


    void plot_equation(        
        char            *   equation,
        char            *   title
    )
    {
        char const *    command    = (handle->nplots > 0) ? "replot" : "plot";
        title                  = (title == NULL)      ? "(none)" : title;

        cmd("%s %s title \"%s\" with %s", command, equation, title, handle->pstyle) ;
        handle->nplots++ ;
        return ;
    }


    int write_x_csv(
        char const * fileName,
        T const * d,
        int n,
        char const * title)
    {
        int     i;
        FILE*   fileHandle;

        if (fileName==NULL || d==NULL || (n<1))
        {
            return -1;
        }

        fileHandle = fopen(fileName, "w");

        if (fileHandle == NULL)
        {
            return -1;
        }

        // Write Comment.
        if (title != NULL)
        {
            fprintf(fileHandle, "# %s\n", title) ;
        }

        /* Write data to this file  */
        for (i=0 ; i<n; i++)
        {
            fprintf(fileHandle, "%d, %.18e\n", i, d[i]) ;
        }

        fclose(fileHandle) ;

        return 0;
    }

    int write_xy_csv(
        char const *        fileName,
        T const    *   x,
        T const    *   y,
        int                 n,
        char const      *   title)
    {
        int     i ;
        FILE*   fileHandle;

        if (fileName==NULL || x==NULL || y==NULL || (n<1))
        {
            return -1;
        }

        fileHandle = fopen(fileName, "w");

        if (fileHandle == NULL)
        {
            return -1;
        }

        // Write Comment.
        if (title != NULL)
        {
            fprintf(fileHandle, "# %s\n", title) ;
        }

        /* Write data to this file  */
        for (i=0 ; i<n; i++)
        {
            fprintf(fileHandle, "%.18e, %.18e\n", x[i], y[i]) ;
        }

        fclose(fileHandle) ;

        return 0;
    }

    int write_multi_csv(
        char const *        fileName,
        T const    **  xListPtr,
        int                 n,
        int                 numColumns,
        char const      *   title)
    {
        int     i;
        int     j;
        FILE*   fileHandle;

        if (fileName==NULL || xListPtr==NULL || (n<1) || numColumns <1)
        {
            return -1;
        }

        for (j=0;j<numColumns;j++)
        {
            if (xListPtr[j] == NULL)
            {
                return -1;
            }
        }

        fileHandle = fopen(fileName, "w");

        if (fileHandle == NULL)
        {
            return -1;
        }

        // Write Comment.
        if (title != NULL)
        {
            fprintf(fileHandle, "# %s\n", title) ;
        }

        /* Write data to this file  */
        for (i=0 ; i<n; i++)
        {
            fprintf(fileHandle, "%d, %.18e", i, xListPtr[0][i]) ;
            for (j=1;j<numColumns;j++)
            {
                fprintf(fileHandle, ", %.18e", xListPtr[j][i]) ;
            }
            fprintf(fileHandle, "\n");
        }

        fclose(fileHandle) ;

        return 0;
    }

    char const * tmpfile()
    {
        static char const * tmp_filename_template = "tmpdatafile_XXXXXX";
        char *              tmp_filename = NULL;
        int                 tmp_filelen = strlen(tmp_filename_template);

    #ifndef WIN32
        int                 unx_fd;
    #endif // #ifndef WIN32

        assert(handle->tmp_filename_tbl[handle->ntmp] == NULL);

        /* Open one more temporary file? */
        if (handle->ntmp == GP_MAX_TMP_FILES - 1) {
            fprintf(stderr,
                    "maximum # of temporary files reached (%d): cannot open more",
                    GP_MAX_TMP_FILES) ;
            return NULL;
        }

        tmp_filename = (char*) malloc(tmp_filelen+1);
        if (tmp_filename == NULL)
        {
            return NULL;
        }
        strcpy(tmp_filename, tmp_filename_template);

    #ifdef WIN32
        if (_mktemp(tmp_filename) == NULL)
        {
            return NULL;
        }
    #else // #ifdef WIN32
        unx_fd = mkstemp(tmp_filename);
        if (unx_fd == -1)
        {
            return NULL;
        }
        close(unx_fd);

    #endif // #ifdef WIN32

        handle->tmp_filename_tbl[handle->ntmp] = tmp_filename;
        handle->ntmp ++;
        return tmp_filename;
    }

    void plot_atmpfile(char const* tmp_filename, char const* title)
    {
        char const *    command    = (handle->nplots > 0) ? "replot" : "plot";
        if( title ){
            cmd("%s \"%s\" title \"%s\" with %s", command, tmp_filename,title, handle->pstyle) ;
        }else{
            cmd("%s \"%s\" notitle with %s", command, tmp_filename,
                    handle->pstyle) ;
        }
        handle->nplots++ ;
        return ;
    }    
};