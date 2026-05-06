using namespace std;
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char **argv)
{

  std::string dir = ".";
  char *versionstr = getenv("USER");
  if (versionstr)
  {
    std::cout << dir + "/VER" + versionstr;
    ofstream version_file(dir + "/VER" + versionstr, fstream::trunc);
    // 	version_file<<versionstr;
    // 	version_file.close();
  }
  else
  {
    std::cout << "HBT_VERSION not specified" << std::endl;
  }

  return 0;
}
