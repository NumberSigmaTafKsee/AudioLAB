%module filesystem
%{
#include <filesystem>    
#include <string>
#include <iostream>

using namespace std::filesystem;
namespace jiminy {}
using namespace jiminy;
%}
%include "std_string.i"
%include "std_vector.i"

%template(string_vector) std::vector<std::string>;

namespace std::filesystem
{
    class path
    {
    public:

        path() noexcept;
        path(const path &p);
        path(const std::string&);

        ~path();

        path& operator = (const path& p);
        path& operator = (const std::string& p);
        //path& operator / (const std::string & s);
        path& append(const std::string &s);
        path& concat(const std::string &s);

        void clear() noexcept;
        path& make_preferred();
        path& remove_filename();
        path& replace_filename(const path& r);
        path& replace_extension(const path& r = path());
        void swap(path & o) noexcept;
        const char* c_str() const noexcept;
        std::string string() const;
        std::wstring wstring() const;
        std::u16string u16string() const;
        std::u32string u32string() const;
        std::string u8string() const;
        //std::u8string u8string() const;
        std::string generic_string() const;
        std::wstring generic_wstring() const;
        std::u16string generic_u16string() const;
        std::u32string generic_u32string() const;
        std::string generic_u8string() const;
        int compare(const path & p) const noexcept;
        path lexically_normal() const;
        path lexically_relative(const path&) const;
        path lexically_proximate(const path&) const;
        path root_name() const;
        path root_directory() const;
        path root_path() const;
        path relative_path() const;
        path parent_path() const;
        path filename() const;
        path stem() const;
        path extension() const;
        bool empty() const;
        bool has_root_path() const;
        bool has_root_name() const;
        bool has_root_directory() const;
        bool has_relative_path() const;
        bool has_parent_path() const;
        bool has_filename() const;
        bool has_stem() const;
        bool has_extension() const;
        bool is_absolute() const;
        bool is_relative() const;

        %extend {
            std::string __str__() { return $self->string(); }
        }
    };

    enum class file_type {
        none,
        not_found,
        regular,
        directory,
        symlink,
        block,
        character,
        fifo,
        socket,
        unknown,
        /* implementation-defined */
    };
    enum class perm_options {        
        replace,
        add,
        remove,
        nofollow
    };
    enum class copy_options {
        none,
        skip_existing,
        overwrite_existing,
        update_existing,
        recursive,
        copy_symlinks,
        skip_symlinks,
        directories_only,
        create_symlinks,
        create_hard_links
    };

    enum class directory_options {
        none,
        follow_directory_symlink,
        skip_permission_denied
    };
		


    struct space_info {
        uintmax_t capacity;
        uintmax_t free;
        uintmax_t available;
    };

    class file_status
    {
    public:
        file_status();
        file_status(std::filesystem::file_status s);

        file_status& operator = (const file_status & s);

        int type() const;
        int permissions() const;
    };

    class directory_entry
    {
    private:
        directory_entry();
        directory_entry(const std::filesystem::directory_entry & e);

        directory_entry& operator=(const directory_entry& e);
        std::string path();
        path to_path();

        void replace_filename(const path&);
        void assign(const path&);
        bool exists() const;
        void refresh();

        bool is_block_file();
        bool is_character_file();
        bool is_directory();
        bool is_fifo();
        bool is_other();
        bool is_regular_file();
        bool is_socket();
        bool is_symlink();
      
        uintmax_t hard_link_count();
        
        //std::filesystem::file_time_type last_write_time() const;
        //std::filesystem::file_time_type last_write_time(std::error_code & ) const;

        file_status status() const;
        
        int type() const;
        int permissions() const;
    };
}

%inline 
{
    std::vector<std::string> dir(const std::string &p)
    {                    
        std::vector<std::string> r;   
        // it crashes here dont know why
        for( auto const& dir_entry: std::filesystem::directory_iterator(p))
        {                                        
            r.push_back(dir_entry.path().string());
        }    
        return r;
    }
    std::vector<std::string> rdir(const std::string &p)
    {                    
        std::vector<std::string> r;                
        for(const auto & dir_entry: std::filesystem::recursive_directory_iterator(p))
        {                
            r.push_back(dir_entry.path().string());
        }
        return r;
    }

    namespace jiminy
    {
    static std::string absolute(const std::string& p) { return std::filesystem::absolute(p).string(); }
    static std::string canonical(const std::string& p) { return std::filesystem::canonical(p).string(); }
    static std::string weakly_canonical(const std::string& p) { return std::filesystem::weakly_canonical(p).string(); }
    static std::string relative(const std::string& p) { return std::filesystem::relative(p).string(); }
    static std::string proximate(const std::string & p) { return std::filesystem::proximate(p).string(); }
    
    static void copy(const std::string & from, const std::string & to) { std::filesystem::copy(from, to); }
    static void copy_file(const std::string & from, const std::string & to) { std::filesystem::copy_file(from, to); }
    static void copy_symlink(const std::string & from, const std::string & to) { std::filesystem::copy_symlink(from, to); }
    
    static bool create_directory(const std::string & path) { return std::filesystem::create_directory(path); }
    static bool create_directories(const std::string & path) { return std::filesystem::create_directories(path); }

    static void create_hard_link(const std::string & target, const std::string & link) { std::filesystem::create_hard_link(target,link); }
    static void create_symlink(const std::string & target, const std::string & link) { std::filesystem::create_symlink(target,link); }
    static void create_directory_symlink(const std::string & target, const std::string & link) { std::filesystem::create_directory_symlink(target,link); }

    static std::filesystem::path current_path() { return (std::filesystem::current_path()); }
    static bool exists(const std::string & p) { return std::filesystem::exists(p); }
    static bool equivalent(const std::string & p1, const std::string & p2) { return std::filesystem::equivalent(p1,p2); }

    static uintmax_t file_size(const std::string & p) { return std::filesystem::file_size(p); }
    static uintmax_t hard_link_count(const std::string & p) { return std::filesystem::hard_link_count(p); }
    
    static std::string read_symlink(const std::string & p) { return std::filesystem::read_symlink(p).string(); }

    static bool remove(const std::string & p) { return std::filesystem::remove(p); }
    static uintmax_t remove_all(const std::string & p) { return std::filesystem::remove_all(p); }        
    static void rename(const std::string & old, const std::string &newf) { std::filesystem::rename(old,newf); }

    static void resize_file(const std::string &p, uintmax_t size) { std::filesystem::resize_file(p, size); }
    static std::string temp_directory_path() { return std::filesystem::temp_directory_path().string(); }

    static bool is_block_file(const std::string &p) { return std::filesystem::is_block_file(p); }
    static bool is_character_file(const std::string &p) { return std::filesystem::is_character_file(p); }
    static bool is_directory(const std::string &p) { return std::filesystem::is_directory(p); }
    static bool is_empty(const std::string &p) { return std::filesystem::is_empty(p); }
    static bool is_fifo(const std::string &p) { return std::filesystem::is_fifo(p); }
    static bool is_other(const std::string &p) { return std::filesystem::is_other(p); }
    static bool is_regular_file(const std::string &p) { return std::filesystem::is_regular_file(p); }
    static bool is_socket(const std::string &p) { return std::filesystem::is_socket(p); }
    static bool is_symlink(const std::string &p) { return std::filesystem::is_symlink(p); }
    }
}