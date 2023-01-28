%module tinydir
%{
#include "tinydir.h"
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <libgen.h>
#include <sys/stat.h>
#include <stddef.h>

%}

typedef char _tinydir_char_t;

typedef struct tinydir_file
{
	_tinydir_char_t path[_TINYDIR_PATH_MAX];
	_tinydir_char_t name[_TINYDIR_FILENAME_MAX];
	_tinydir_char_t *extension;
	int is_dir;
	int is_reg;
	//struct stat _s;
} 
tinydir_file;

typedef struct tinydir_dir
{
	_tinydir_char_t path[_TINYDIR_PATH_MAX];
	int has_next;
	size_t n_files;
	tinydir_file *_files;
	_TINYDIR_DIR *_d;
	struct _tinydir_dirent *_e;
	//struct _tinydir_dirent *_ep;
} tinydir_dir;


int tinydir_open(tinydir_dir *dir, const _tinydir_char_t *path);
int tinydir_open_sorted(tinydir_dir *dir, const _tinydir_char_t *path);
void tinydir_close(tinydir_dir *dir);
int tinydir_next(tinydir_dir *dir);
int tinydir_readfile(const tinydir_dir *dir, tinydir_file *file);
int tinydir_readfile_n(const tinydir_dir *dir, tinydir_file *file, size_t i);
int tinydir_open_subdir_n(tinydir_dir *dir, size_t i);
int tinydir_file_open(tinydir_file *file, const _tinydir_char_t *path);
void _tinydir_get_ext(tinydir_file *file);
int _tinydir_file_cmp(const void *a, const void *b);




