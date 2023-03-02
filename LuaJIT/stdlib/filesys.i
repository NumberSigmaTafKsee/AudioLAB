%module fs 
%{ 
#include "StdFileSystem.h"
using namespace Std;
%}

%include "stdint.i"
%include "std_common.i"
%include "std_string.i"
%include "std_vector.i"

%include "StdFileSystem.h"

%template(string_vector) std::vector<std::string>;
%template(dir_vector) std::vector<DirEntry>;
%template(path_vector) std::vector<Path>;

%constant int copy_none = std::filesystem::copy_options::none;
%constant int copy_skip_existing = std::filesystem::copy_options::skip_existing;
%constant int copy_overwrite_existing = std::filesystem::copy_options::overwrite_existing;
%constant int copy_recursive = std::filesystem::copy_options::recursive;
%constant int copy_copy_symlinks = std::filesystem::copy_options::copy_symlinks;
%constant int copy_skip_symlinks = std::filesystem::copy_options::skip_symlinks;
%constant int copy_directories_only = std::filesystem::copy_options::directories_only;
%constant int copy_create_symlinks = std::filesystem::copy_options::create_symlinks;
%constant int copy_create_hard_links = std::filesystem::copy_options::create_hard_links;

%constant int file_type_none = std::filesystem::file_type::none;
%constant int file_type_not_found = std::filesystem::file_type::not_found;
%constant int file_type_regular = std::filesystem::file_type::regular;
%constant int file_type_directory = std::filesystem::file_type::directory;
%constant int file_type_symlink = std::filesystem::file_type::symlink;
%constant int file_type_block = std::filesystem::file_type::block;
%constant int file_type_character = std::filesystem::file_type::character;
%constant int file_type_fifo = std::filesystem::file_type::fifo;
%constant int file_type_socket = std::filesystem::file_type::socket;
%constant int file_type_unknown = std::filesystem::file_type::unknown;

%constant int perm_none = std::filesystem::perms::none;
%constant int perm_owner_read = std::filesystem::perms::owner_read;
%constant int perm_owner_write = std::filesystem::perms::owner_write;
%constant int perm_owner_exec = std::filesystem::perms::owner_exec;
%constant int perm_owner_all = std::filesystem::perms::owner_all;
%constant int perm_group_all = std::filesystem::perms::group_all;
%constant int perm_group_read = std::filesystem::perms::group_read;
%constant int perm_group_write = std::filesystem::perms::group_write;
%constant int perm_group_exec = std::filesystem::perms::group_exec;
%constant int perm_others_read = std::filesystem::perms::others_read;
%constant int perm_others_write = std::filesystem::perms::others_write;
%constant int perm_others_exec = std::filesystem::perms::others_exec;
%constant int perm_others_all = std::filesystem::perms::others_all;
%constant int perm_set_uid = std::filesystem::perms::set_uid;
%constant int perm_set_gid = std::filesystem::perms::set_gid;
%constant int perm_stick_bit = std::filesystem::perms::sticky_bit;
%constant int perm_mask = std::filesystem::perms::mask;
%constant int perm_unknown = 0xFFFF;
