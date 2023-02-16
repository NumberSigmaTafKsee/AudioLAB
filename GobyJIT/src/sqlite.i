%module sqlite 
%{
#include <sqlite3pp.h>
#include <cassert>
%}

%include "stdint.i"
%include "lua_fnptr.i"


%import <sqlite3pp.h>

class null_type {};
enum copy_semantic { copy, nocopy };

%inline 
%{
    struct Database;
    // dont use with threads which is not a problem normally with Lua.
    static Database *current_database = NULL;
    static busy_handler_cb(void * db);
    static void backup_handler_cb(void * db, int a, int b, int c);
    static int commit_handler_cb(void * db);
    static void update_handler_cb(void* p, 
                                    int opcode, 
                                    char const* dbname, 
                                    char const* tablename, 
                                    long long int rowid);
    static int authorize_handler_cb(void * db, int evcode, char const* p1, char const* p2, char const* dbname, char const* tvname);

    struct DatabaseError
    {
        sqlite3pp::database_error error;

        DatabaseError(sqlite3pp::database_error & err) {
            error = err;
        }
        
        std::string what() { return error.what(); }    
    };

    struct Database 
    {
        sqlite3pp::database *p_db;
        
        SWIGLUA_FN busy_handler;
        SWIGLUA_FN commit_handler;
        SWIGLUA_FN rollback_handler;

        Database() p_db(NULL) {} 

        Database(const char * dbname = NULL, int flags = SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, const char* vfs = NULL)
        {
            p_db = new database(dbname,flags,vfs);
            assert(p_db != NULL);
        }

        Database(Database& dbase)
        {
            p_db = new database(*dbase.p_db);
            assert(p_db != NULL);
        }

        ~Database() 
        {
            if(p_db) delete p_db;
        }



        int connect(char const* dbname, int flags, const char* vfs = NULL)
        {
            return p_db->connect(dbname,flags,vfs);
        }
        int disconnect()
        {
            return p_db->disconnect();
        }

        int attach(char const* dbname, char const* name)
        {
            return p_db->attach(dbname,name);
        }
        int detach(char const* name)
        {
            return p_db->detach(name);
        }

        
        int backup(database& destdb)
        {
            assert(backup_handler.L != NULL);
            // if using thread need to lock this.
            current_database = this;
            p_db->backup(destdb, backup_handler_cb);
        }
        int backup(char const* dbname, Database& destdb, char const* destdbname, backup_handler h, int step_page = 5)
        {
            assert(backup_handler.L != NULL);
            current_database = this;
            p_db->backup(dbname,*destdb.p_db, destdbname, backup_handler_cb,step_page);
        }
        
        long long int last_insert_rowid() const
        {
            return p_db->last_insert_rowid();
        }

        int enable_foreign_keys(bool enable = true)
        {
            return p_db->enable_foreign_keys(enable);
        }
    
        int enable_triggers(bool enable = true)
        {
            return p_db-=>enable_triggers(enable);
        }
    
        int enable_extended_result_codes(bool enable = true)
        {
            return p_db-> enable_extended_result_codes(enable);
        }

        int changes() const
        {
            return p_db->changes();
        }

        int error_code() const
        {
            return p_db->error_code();
        }
    
        int extended_error_code() const
        {
            return p_db->extended_error_code();
        }

        char const* error_msg() const
        {
            return p_db->error_msg();
        }

        int execute(char const* sql)
        {
            return p_db->execute(sql);
        }

       // int executef(char const* sql, ...);

        
        int set_busy_timeout(int ms)
        {
            return p_db->set_busy_timeout(ms);        
        }

        void set_backup_handle(SWIGLUA_FN h)        
        {
            backup_handler = h;
        }
        void set_busy_handle(SWIGLUA_FN h)        
        {
            busy_handler = h;
        }
        void set_commit_handler(SWIGLUA_FN h)
        {
            commit_handler = h;
        }
        void set_rollback_handler(SWIGLUA_FN h)
        {
            rollback_handler = h;
        }
        void set_update_handler(SWIGLUA_FN h h)
        {
            update_handler = h;
        }
        void set_authorize_handler(SWIGLUA_FN h)
        {
            authorize_handler = h;
        }

    };

    struct Statement
    {
        sqlite3pp::statement * p_st;
        Database db;

        Statement(Database & dbase, char const * stmt = NULL)
        {
            db = dbase;
            p_st = new sqlite3pp::statment(dbase.p_db,stmt);
            assert(p_st != NULL);
        }
        ~Statement() 
        {
            if(p_st) delete p_st;
        }

        int prepare(char const* stmt)
        {
            return p_st->prepare(stmt);
        }

        int finish()
        {
            return p_st->finish();
        }

        int bind(int idx, int value)
        {
            return p_st->bind(idx,value);
        }
        
        int bind(int idx, double value)
        {
            return p_st->bind(idx,value);
        }
        int bind(int idx, long long int value)
        {
            return p_st->bind(idx,value);
        }
        int bind(int idx, char const* value, copy_semantic fcopy) {
            return p_st->bind(idx,value,fcopy);
        }
        int bind(int idx, void const* value, int n, copy_semantic fcopy) {
            return p_st->bind(idx,value,n,fcopy);
        }
        int bind(int idx, std::string const& value, copy_semantic fcopy) {
            return p_st->bin(idx,value,fcopy);
        }
    
        int bind(int idx)
        {
            return p_st->bind(idx);
        }        
        int bind(int idx)
        {
            null_type null;
            return p_st->bind(idx,null);
        }
        int bind(char const* name, int value)
        {
            return p_st->bind(name,value);
        }
        int bind(char const* name, double value)
        {
            return p_st->bind(name,value);
        }
        int bind(char const* name, long long int value)
        {
            return p_st->bind(name,value);
        }
        int bind(char const* name, char const* value, copy_semantic fcopy) {
            return p_st->bind(name,value,fcopy);
        }
        int bind(char const* name, void const* value, int n, copy_semantic fcopy)
        {
            return p_st->bind(name,value,fcopy);
        }
        int bind(char const* name, std::string const& value, copy_semantic fcopy) {
            return p_st->bind(name,value,fcopy);
        }
        
        int bind(char const* name)
        {
            return p_st->bind(name);
        }
        int bind(char const* name) {
            null_type null;
            return p_st->bind(name,null);
        }

        int step() { return p_st->step(); }
        int reset() { return p_st->reset(); }
    };

    
    struct Command 
    {
        sqlite3pp::command *cmd;
        Database db;
        sqlite3pp::command::bindstream binder;

        Command(Database &db, char const * stmt = nullptr)
        {
            this->db = db;
            cmd = new sqlite3pp::command(db->p_db,stmt);
            assert(cmd != NULL);
        }
        ~Command() {
            if(cmd) delete cmd;
        }

        void binder(int idx=1) {
            binder = cmd->binder(idx);
        }
        int execute() {
            cmd->execute();
        }
        int execute_all() {
            cmd->execute_all();
        }

        
        void insert(double value) {
            *binder << value;
        }                
        void insert(int64_t value) {
            *binder << value;
        }        
        void insert(const char * value) {
            *binder << value;
        }
        void insert(const std::string & value) {
            *binder << value;
        }
    };
    typedef struct Query 
    {
        sqlite3pp::query * query;
        sqlite3pp::query_iterator it;

        Query(Statement & stmt) {
            query = new sqlite3pp::query(stmt.p_st);
            assert(query != NULL);
            it = sqlite3pp::query_iterator(query);
        }
        ~Query() {
            if(query) delete query;
        }


        int data_count() const { return query->data_count(); }
        int column_type() const { return query->column_type(); }
        int column_bytes(int idx) const { return query->column_bytes(idx); }

        
        double get_double(int idx) {
            return query->get(idx,(double)0);            
        }
        int64_t get_int64(int idx) {
            return query->get(idx,(int64_t)0);            
        }
        char const * get_cstr(int idx) {
            return query->get(idx, (char const *)0);
        }
        void const * get_ptr(int idx) {
            return query->get(idx, (void const*)0);
        }
        null_type get_null(int idx) {
            null_type null;
            return query->get(idx,null);
        }
                
    };

    struct Transaction
    {
        sqlite3pp::transaction *tr;

        Transaction(Database & db, bool fcommit=false, bool freserve = false) {
            tr = new sqlite3pp::transaction(db.p_db,fcommit,freserve);
        }
        ~Transaction() {
            if(tr) delete tr;
        }
        int commint() { return tr->commit(); }
        int rollback() { return tr->rollback(); }        
    };
    

    static int authorize_handler_cb(void * db, int evcode, char const* p1, char const* p2, char const* dbname, char const* tvname)
    {
        Database * p = (Database*)db;   
        assert(p != NULL);
        assert(p->busy_handler.L != NULL);
        SWIGLUA_FN_GET(p->busy_handler);
        lua_pushnumber(evcode);
        lua_pushstring(p1);
        lua_pushstring(p2);
        lua_pushstring(dbname);
        lua_pushstring(tvname);
        assert(lua_pcall(p->busyhandler.L,1,5,1) == 0);        
        return lua_tonumber(p->busyhandler.L,-1);
    }

    static int busy_handler_cb(void * db, int cnt)
    {
        Database * p = (Database*)db;   
        assert(p != NULL);
        assert(p->busy_handler.L != NULL);
        SWIGLUA_FN_GET(p->busy_handler);
        lua_pushnumber(cnt);
        assert(lua_pcall(p->busyhandler.L,1,1,0) == 0);        
        return lua_tonumber(p->busyhandler.L,-1);
    }

    static void backup_handler_cb(void * db, int a, int b, int c)
    {
        // db will be ptr to database which is already in current_database.
        Database * p = current_database;
        assert(p != NULL);
        assert(p->busy_handler.L != NULL);
        SWIGLUA_FN_GET(p->busy_handler);
        lua_pushnumber(a);
        lua_pushnumber(b);
        lua_pushnumber(c);
        assert(lua_pcall(p->busyhandler.L,3,0,0) == 0);                
    }

    static int commit_handler_cb(void * db)
    {
        Database * p = (Database*)db;   
        assert(p != NULL);
        assert(p->busy_handler.L != NULL);
        SWIGLUA_FN_GET(p->busy_handler);        
        assert(lua_pcall(p->busyhandler.L,0,1,0) == 0);        
        return lua_tonumber(p->busyhandler.L,-1);
    }

    static void update_handler_cb(void* p, 
                                    int opcode, 
                                    char const* dbname, 
                                    char const* tablename, 
                                    long long int rowid)
    {
        Database * p = (Database*)db;   
        assert(p != NULL);
        assert(p->busy_handler.L != NULL);
        SWIGLUA_FN_GET(p->busy_handler);  
        lua_pushnumber(opcode);
        lua_pushstring(dbname);
        lua_pushstring(tablename);
        lua_pushstring(rowid);
        assert(lua_pcall(p->busyhandler.L,0,4,0) == 0);        
        return lua_tonumber(p->busyhandler.L,-1);
    }
                            
%}
