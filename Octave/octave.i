////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 1996-2022 The Octave Project Developers
//
// See the file COPYRIGHT.md in the top-level directory of this
// distribution or <https://octave.org/copyright/>.
//
// This file is part of Octave.
//
// Octave is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Octave is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Octave; see the file COPYING.  If not, see
// <https://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////

%module octopus
%{
#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>
#include <octave/interpreter.h>
using namespace octave;
%}

//%include "Octave/std_math.i"
%include "std_vector.i"
%include "std_string.i"

// if you dont have a problem uncomment this
//%include "Octave/stdint.i"
// and comment out these
typedef signed char int8_t;
typedef unsigned char uint8_t;
typedef signed short int16_t;
typedef unsigned short uint16_t;
typedef signed int int32_t;
typedef unsigned int uint32_t;
typedef signed long int64_t;
typedef unsigned long uint64_t;

%include <octave/oct.h>
#undef octave_octave_config_h
#undef OCTAVE_AUTOCONFIG_H_INCLUDED
%include <octave/octave-config.h>
%include <octave/oct-inttypes-fwd.h>
%include <octave/oct-cmplx.h>

/*
/snap/octave/current/include/octave-7.1.0/octave/oct-atomic.h                    /snap/octave/current/include/octave-7.1.0/octave/oct-locbuf.h
/snap/octave/current/include/octave-7.1.0/octave/octave-build-info.h             /snap/octave/current/include/octave-7.1.0/octave/oct-lvalue.h
/snap/octave/current/include/octave-7.1.0/octave/octave-config.h                 /snap/octave/current/include/octave-7.1.0/octave/oct-map.h
/snap/octave/current/include/octave-7.1.0/octave/octave-default-image.h          /snap/octave/current/include/octave-7.1.0/octave/oct-mutex.h
/snap/octave/current/include/octave-7.1.0/octave/octave.h                        /snap/octave/current/include/octave-7.1.0/octave/oct-norm.h
/snap/octave/current/include/octave-7.1.0/octave/octave-preserve-stream-state.h  /snap/octave/current/include/octave-7.1.0/octave/oct-password.h
/snap/octave/current/include/octave-7.1.0/octave/oct-base64.h                    /snap/octave/current/include/octave-7.1.0/octave/oct-prcstrm.h
/snap/octave/current/include/octave-7.1.0/octave/oct-binmap.h                    /snap/octave/current/include/octave-7.1.0/octave/oct-procbuf.h
/snap/octave/current/include/octave-7.1.0/octave/oct-cmplx.h                     /snap/octave/current/include/octave-7.1.0/octave/oct-process.h
/snap/octave/current/include/octave-7.1.0/octave/oct-convn.h                     /snap/octave/current/include/octave-7.1.0/octave/oct-rand.h
/snap/octave/current/include/octave-7.1.0/octave/oct-env.h                       /snap/octave/current/include/octave-7.1.0/octave/oct-refcount.h
/snap/octave/current/include/octave-7.1.0/octave/oct-errno.h                     /snap/octave/current/include/octave-7.1.0/octave/oct-rl-edit.h
/snap/octave/current/include/octave-7.1.0/octave/oct-fftw.h                      /snap/octave/current/include/octave-7.1.0/octave/oct-rl-hist.h
/snap/octave/current/include/octave-7.1.0/octave/oct-fstrm.h                     /snap/octave/current/include/octave-7.1.0/octave/oct-shlib.h
/snap/octave/current/include/octave-7.1.0/octave/oct-glob.h                      /snap/octave/current/include/octave-7.1.0/octave/oct-sort.h
/snap/octave/current/include/octave-7.1.0/octave/oct-group.h                     /snap/octave/current/include/octave-7.1.0/octave/oct-spparms.h
/snap/octave/current/include/octave-7.1.0/octave/oct.h                           /snap/octave/current/include/octave-7.1.0/octave/oct-stdstrm.h
/snap/octave/current/include/octave-7.1.0/octave/oct-handle.h                    /snap/octave/current/include/octave-7.1.0/octave/oct-stream.h
/snap/octave/current/include/octave-7.1.0/octave/oct-hdf5-types.h                /snap/octave/current/include/octave-7.1.0/octave/oct-string.h
/snap/octave/current/include/octave-7.1.0/octave/oct-hist.h                      /snap/octave/current/include/octave-7.1.0/octave/oct-strstrm.h
/snap/octave/current/include/octave-7.1.0/octave/oct-inttypes-fwd.h              /snap/octave/current/include/octave-7.1.0/octave/oct-syscalls.h
/snap/octave/current/include/octave-7.1.0/octave/oct-inttypes.h                  /snap/octave/current/include/octave-7.1.0/octave/oct-time.h
/snap/octave/current/include/octave-7.1.0/octave/oct-iostrm.h                    /snap/octave/current/include/octave-7.1.0/octave/oct-uname.h
*/

%include "Octave/dim_vector.i"
%include "Octave/array.i"
%include "Octave/marray.i"
%include "Octave/cNDArray.i"
%include "Octave/cRowVector.i"
%include "Octave/cColVector.i"
%include "Octave/cMatrix.i"
%include "Octave/dNDArray.i"
%include "Octave/dRowVector.i"
%include "Octave/dColVector.i"
%include "Octave/dMatrix.i"
%include "Octave/fNDArray.i"
%include "Octave/fRowVector.i"
%include "Octave/fColVector.i"
%include "Octave/fMatrix.i"
%include "Octave/fCNDArray.i"
%include "Octave/fCRowVector.i"
%include "Octave/fCColVector.i"
%include "Octave/fCMatrix.i"
%include "Octave/octave_value.i"
%include "Octave/octave_value_list.i"
//%include "Octave/octave_cell.i"

namespace octave
{
    ///////////////////////////////
    class interpreter
    {
    public:

        // Create an interpreter object and perform basic initialization.

        interpreter ();

        // No copying, at least not yet...
        interpreter (const interpreter&) = delete;
        interpreter& operator = (const interpreter&) = delete;

        // Clean up the interpreter object.
        ~interpreter (void);

        void intern_nargin (octave_idx_type nargs);

        // If creating an embedded interpreter, you may inhibit reading
        // the command history file by calling initialize_history with
        // read_history_file = false prior to calling initialize.

        void initialize_history (bool read_history_file = false);

        // If creating an embedded interpreter, you may inhibit setting
        // the default compiled-in path by calling initialize_load_path
        // with set_initial_path = false prior calling initialize.  After
        // that, you can add directories to the load path to set up a
        // custom path.

        void initialize_load_path (bool set_initial_path = true);

        // Load command line history, set the load path.

        void initialize (void);

        // Note: GET_LINE_AND_EVAL is only used by new experimental terminal
        // widget.

        void get_line_and_eval (void);

        // Parse a line of input.  If input ends at a complete statement
        // boundary, execute the resulting parse tree.  Useful to handle
        // parsing user input when running in server mode.

        void parse_and_execute (const std::string& input, bool& incomplete_parse);

        // Initialize the interpreter (if not already done by an explicit
        // call to initialize), execute startup files, --eval option code,
        // script files, and/or interactive commands.

        int execute (void);

        bool server_mode (void) const;

        bool interactive (void) const;

        void interactive (bool arg);
        void read_site_files (bool flag);
        void read_init_files (bool flag);
        void verbose (bool flag);
        void traditional (bool flag);
        bool traditional (void) const;
        void inhibit_startup_message (bool flag);
        bool in_top_level_repl (void) const;
        bool initialized (void) const;
        void interrupt_all_in_process_group (bool b);
        bool interrupt_all_in_process_group (void) const;
        /*
        application *get_app_context (void);
        display_info& get_display_info (void);
        environment& get_environment (void);
        settings& get_settings (void);
        error_system& get_error_system (void);
        help_system& get_help_system (void);
        input_system& get_input_system (void);
        output_system& get_output_system (void);
        history_system& get_history_system (void);
        dynamic_loader& get_dynamic_loader (void);
        load_path& get_load_path (void);
        load_save_system& get_load_save_system (void);
        type_info& get_type_info (void);
        symbol_table& get_symbol_table (void);
        tree_evaluator& get_evaluator (void);
        
        symbol_scope get_top_scope (void) const;
        symbol_scope get_current_scope (void) const;
        symbol_scope require_current_scope (const std::string& who) const;

        profiler& get_profiler (void);

        stream_list& get_stream_list (void);

        child_list& get_child_list (void);
        url_handle_manager& get_url_handle_manager (void);

        cdef_manager& get_cdef_manager (void);
        gtk_manager& get_gtk_manager (void);
        event_manager& get_event_manager (void);
        gh_manager& get_gh_manager (void);
        */
        // Any Octave code that needs to change the current directory should
        // call this function instead of calling the system chdir function
        // directly so that the  load-path and GUI may be notified of the
        // change.

        int chdir (const std::string& dir);

        void mlock (bool skip_first = false) const;
        void munlock (bool skip_first = false) const;
        bool mislocked (bool skip_first = false) const;

        // NOTE: since we have a version that accepts a bool argument, we
        // can't rely on automatic conversion from char* to std::string.
        //void munlock (const char *nm);
        void munlock (const std::string& nm);

        //bool mislocked (const char *nm);
        bool mislocked (const std::string& nm);

        std::string mfilename (const std::string& opt = "") const;

        octave_value_list eval_string (const std::string& eval_str, bool silent,
                                    int& parse_status, int nargout);

        octave_value eval_string (const std::string& eval_str, bool silent,
                                int& parse_status);

        octave_value_list eval_string (const octave_value& arg, bool silent,
                                    int& parse_status, int nargout);

        octave_value_list eval (const std::string& try_code, int nargout);

        octave_value_list eval (const std::string& try_code,
                                const std::string& catch_code, int nargout);

        octave_value_list evalin (const std::string& context,
                                const std::string& try_code, int nargout);

        octave_value_list evalin (const std::string& context,
                                const std::string& try_code,
                                const std::string& catch_code, int nargout);

        /*
        octave_value_list
        feval (const char *name,
            const octave_value_list& args = octave_value_list (),
            int nargout = 0);
        */
        octave_value_list
        feval (const std::string& name,
            const octave_value_list& args = octave_value_list (),
            int nargout = 0);

        octave_value_list
        feval (octave_function *fcn,
            const octave_value_list& args = octave_value_list (),
            int nargout = 0);

        octave_value_list
        feval (const octave_value& f_arg,
            const octave_value_list& args = octave_value_list (),
            int nargout = 0);

        octave_value_list feval (const octave_value_list& args, int nargout = 0);


        octave_value make_function_handle (const std::string& name);
        
        void install_variable (const std::string& name, const octave_value& value,
                            bool global);

        //void set_global_value (const std::string& name, const octave_value& value);

        octave_value global_varval (const std::string& name) const;

        void global_assign (const std::string& name,
                            const octave_value& val = octave_value ());

        octave_value top_level_varval (const std::string& name) const;

        void top_level_assign (const std::string& name,
                            const octave_value& val = octave_value ());

        bool is_variable (const std::string& name) const;

        bool is_local_variable (const std::string& name) const;

        octave_value varval (const std::string& name) const;

        void assign (const std::string& name,
                    const octave_value& val = octave_value ());

        void assignin (const std::string& context, const std::string& varname,
                    const octave_value& val = octave_value ());

        void source_file (const std::string& file_name,
                        const std::string& context = "",
                        bool verbose = false, bool require_file = true);

        bool at_top_level (void) const;

        bool isglobal (const std::string& name) const;

        octave_value find (const std::string& name);

        void clear_all (bool force = false);

        void clear_objects (void);

        void clear_variable (const std::string& name);

        void clear_variable_pattern (const std::string& pattern);

        void clear_variable_regexp (const std::string& pattern);

        void clear_variables (void);

        void clear_global_variable (const std::string& name);

        void clear_global_variable_pattern (const std::string& pattern);

        void clear_global_variable_regexp (const std::string& pattern);

        void clear_global_variables (void);

        void clear_functions (bool force = false);

        void clear_function (const std::string& name);

        void clear_symbol (const std::string& name);

        void clear_function_pattern (const std::string& pat);

        void clear_function_regexp (const std::string& pat);

        void clear_symbol_pattern (const std::string& pat);

        void clear_symbol_regexp (const std::string& pat);

        std::list<std::string> variable_names (void);

        std::list<std::string> top_level_variable_names (void);

        std::list<std::string> global_variable_names (void);

        std::list<std::string> user_function_names (void);

        std::list<std::string> autoloaded_functions (void) const;

        void interrupt (void);

        // Pause interpreter execution at the next available statement and
        // enter the debugger.
        void pause (void);

        // Exit debugger or stop execution and return to the top-level REPL
        // or server loop.
        void stop (void);

        // Add EXPR to the set of expressions that may be evaluated when the
        // debugger stops at a breakpoint.
        void add_debug_watch_expression (const std::string& expr);

        // Remove EXPR from the set of expressions that may be evaluated
        // when the debugger stops at a breakpoint.
        void remove_debug_watch_expression (const std::string& expr);

        // Clear the set of expressions that may be evaluated when the
        // debugger stops at a breakpoint.
        void clear_debug_watch_expressions (void);

        // Return the set of expressions that may be evaluated when the
        // debugger stops at a breakpoint.
        std::set<std::string> debug_watch_expressions (void) const;

        // Resume interpreter execution if paused.
        void resume (void);

        // Provided for convenience.  Will be removed once we eliminate the
        // old terminal widget.
        bool experimental_terminal_widget (void) const;

        //void handle_exception (const execution_exception& ee);

        void recover_from_exception (void);

        void mark_for_deletion (const std::string& file);

        void cleanup_tmp_files (void);

        void quit (int exit_status, bool force = false, bool confirm = true);

        void cancel_quit (bool flag);

        bool executing_finish_script (void) const;

        void add_atexit_fcn (const std::string& fname);

        bool remove_atexit_fcn (const std::string& fname);

        static interpreter * the_interpreter (void) { return m_instance; }

    };

    extern std::string
    get_help_from_file (const std::string& nm, bool& symbol_found,
                        std::string& file);

    extern std::string
    get_help_from_file (const std::string& nm, bool& symbol_found);

    extern octave_value
    load_fcn_from_file (const std::string& file_name,
                        const std::string& dir_name = "",
                        const std::string& dispatch_type = "",
                        const std::string& package_name = "",
                        const std::string& fcn_name = "",
                        bool autoload = false);

    extern void
    source_file (const std::string& file_name,
                const std::string& context = "",
                bool verbose = false, bool require_file = true);

    /*
    extern octave_value_list
    feval (const char *name,
            const octave_value_list& args = octave_value_list (),
            int nargout = 0);
    */
    extern octave_value_list
    feval (const std::string& name,
            const octave_value_list& args = octave_value_list (),
            int nargout = 0);

    extern octave_value_list
    feval (octave_function *fcn,
            const octave_value_list& args = octave_value_list (),
            int nargout = 0);

    extern octave_value_list
    feval (const octave_value& val,
            const octave_value_list& args = octave_value_list (),
            int nargout = 0);

    extern octave_value_list
    feval (const octave_value_list& args, int nargout = 0);

    //extern void
    //cleanup_statement_list (tree_statement_list **lst);
}

//typename Alloc = std::pmr::polymorphic_allocator<T>;
//%template (array_double) Array<double,std::pmr::polymorphic_allocator<T>>;
//%template (array_double) Array<double,std::allocator<T>>;
//%template (marray_double) MArray<double>;
