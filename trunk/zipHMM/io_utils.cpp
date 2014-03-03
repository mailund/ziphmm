#include "io_utils.hpp"

#if not defined WIN32
#include <unistd.h>
#include <pwd.h>
#include <sys/stat.h>
#include <sys/types.h>
#else
#include <windows.h>
#endif


namespace zipHMM {
  namespace {
    void checkTokenOrExit(const std::string &expected, const std::string &found) {
      if(expected != found) {
	std::cerr << "Expected to read token \"" << expected
		  << "\" but got \"" << found << "\"" << std::endl;
	exit(-1);
      }
    }
  }

  bool readToken(std::istream &in, const std::string &expected) {
    std::string found;
    if(!(in >> found))
      return false;
    checkTokenOrExit(expected, found);
    return true;
  }

  void read_token_or_exit(std::istream &in, const std::string &expected) {
    if(!readToken(in, expected)) {
      std::cerr << "Error trying to read token " << expected << std::endl;
      exit(-1);
    }
  }

  bool read_token_or_tell(std::istream &in, const std::string &expected) {
    std::string found;
    if(!(in >> found)) {
      std::cerr << "Error trying to read token " << expected << std::endl;
      exit(-1);
    }
    if(found == expected) {
      return true;
    }

    // else rewind
    in.seekg(in.tellg() - std::streamoff(found.length()), std::ios::beg);
    return false;
  }


  std::string get_user_dir() {
    const char *dir = 0;
#if defined WIN32
    dir = getenv("APPDATA");
    if(dir != 0)
      return dir;

    dir = getenv("USERPROFILE");
    if(dir != 0)
      return dir;

    return std::string(getenv("HOMEDRIVE")) + std::string(getenv("HOMEPATH"));
#else
    dir = getenv("HOME");
    if(dir != 0)
      return dir;

    passwd *pwd = getpwuid(getuid());
    return pwd->pw_dir;
#endif
  }

  std::string get_working_directory() {
    char temp[MAX_PATH_LENGTH];
    return ( getcwd(temp, MAX_PATH_LENGTH) ? std::string( temp ) : std::string("") );
  }

  void mk_dir(const std::string &path) {
#if defined WIN32
    _mkdir(path.c_str());
#else 
    mkdir(path.c_str(), 0777);
#endif
  }

}
