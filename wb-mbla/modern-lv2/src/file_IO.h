// Gunwoo Park's COND interface, modernized to value-owned validated settings.
#ifndef READ_FILE_CONDITION_H
#define READ_FILE_CONDITION_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <map>
#include <stdexcept>
#define TRUE 1
#define FALSE 0
using namespace std;
class COND {
  std::map<std::string, std::string> values_;
public:
  explicit COND(const char* filename);
  const std::string& operator()(const std::string& key) const;
  int cond_print() const;
  void validate() const;
};
// An opened stream can still fail later (quota, disk or close errors).
// Call after writing each owned output, before reporting successful execution.
void finish_output(std::ofstream& output, const std::string& path);
#endif
