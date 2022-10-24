#ifndef PAEAN_UTIL_H
#define PAEAN_UTIL_H

#include <string>
#include <vector>
#include <sstream>
#include <cstring>
#include <cstdint>
#include <unistd.h>
#include <sys/stat.h>
#include <initializer_list>

static inline bool isDigitstr(char *str) 
{
    return strspn(str, "0123456789") == strlen(str);
}

static inline unsigned ceilDiv(unsigned a, unsigned b)
{
    return (a + b - 1) / b;
}

static inline uint32_t nextPow2(uint32_t x)
{
    x--;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x++;
    return x;
}

static inline bool endsWith(std::string const &fullString, std::string const &ending) 
{
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare(fullString.length() - ending.length(), 
                                        ending.length(), ending));
    } else {
        return false;
    }
}

static inline bool eraseSuffixStr(std::string &fullString, std::string const &subString)
{
    if (endsWith(fullString, subString)) {
        fullString = fullString.substr(0, fullString.size() - subString.size());
        return true;
    } else {
        return false;
    }
}

static inline bool isdir(const char *dir)
{
    struct stat buf;
    if (stat(dir, &buf) == 0 && S_ISDIR(buf.st_mode))
        return true;
    return false;
}

static inline void createDirIfNotExists(const char *dir)
{
    if (!isdir(dir)) {
        mkdir(dir, 0777);
    }
}

static inline std::string pathJoin(std::string path1, std::string path2)
{   
    if (endsWith(path1, "/")) 
        return path1 + path2;
    else
        return path1 + "/" + path2;
}

static inline std::string extractFileName(char *file_path)
{
    std::string path = std::string(file_path);
    std::size_t last = path.find_last_of("/");
    if (last != std::string::npos) {
        return path.substr(last+1, path.size()-last-1);
    } else {
        return path;
    }
}

static inline std::string extractFilePrefix(char *file_path)
{
    std::string file = extractFileName(file_path);
    std::size_t last = file.find_last_of(".");
    if (last != std::string::npos) {
        return file.substr(0, last);
    } else {
        return file;
    }
}

static inline std::string getPwd()
{
    char buffer[PATH_MAX];
    char *res = getcwd(buffer, PATH_MAX);
    
    // ignore compiler warning
    (void)res;

    return std::string(buffer); 
}

static inline int catoi(char *s) 
{
    int acum = 0;
    while((*s >= '0')&&(*s <= '9')) {
        acum = acum * 10;
        acum = acum + (*s - 48);
        s++;
    }
    return acum;
}

static inline void trim(std::string& s, std::string& t)
{
    s.erase(s.find_last_not_of(t) + 1);
    s.erase(0, s.find_first_not_of(t));
}

static inline std::string join(std::initializer_list<std::string> strings,
                               char delimiter=';')
{
    std::stringstream ss;
    for (auto it = strings.begin(); it != strings.end(); it++) {
        if (it != (strings.end() - 1))
            ss << *it << delimiter;
        else
            ss << *it;
    }
    return ss.str();
}

static inline std::vector<std::string> split(const std::string &s,
                                             char delimiter=';') {
    std::stringstream ss(s);
    std::string str;
    std::vector<std::string> strings;
    while (std::getline(ss, str, delimiter)) {
        strings.emplace_back(str);
    }
    return strings;
}

#endif //PAEAN_UTIL_H
