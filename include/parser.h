#include <string>
#include <fstream>
#include <algorithm>

#ifndef PARSER_H
#define	PARSER_H
  // Parse a configuration file
  // @param   fileName The name of the file to parse
  // @return  none
  void ParserConfigFile(std::string fileName);
  std::string GetOptionToString(std::string key);
  const char *GetOptionToChar(std::string key);
  int GetOptionToInt(std::string key);
  double GetOptionToDouble(std::string key);
  
  std::string ParseName(std::string value);
  std::string ParseValue(std::string value);
  std::string Trim(std::string str);
  std::string RTrim(std::string str);
  std::string LTrim(std::string str);

#endif	// PARSER_H