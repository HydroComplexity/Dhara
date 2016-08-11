/*
// Copyright (C) 2016, HydroComplexity Group
// All rights reserved.
//
// Distributed Hydrologicc and Regional Analysis (DHARA) Model
// DHARA model is made available as a restricted, non-exclusive, 
// non-transferable license for education and research purpose only, 
// and not for commercial use. See the LICENSE.txt for more details.
//
// Author: levuvietphong@gmail.com (Phong Le)
*/


#include <string.h>
#include <functional>
#include "../include/main.h"
#include "../include/parser.h"


////////////////////////////////////////////////////////
// PARSING CONFIGURATION FILE READ FROM COMMAND LINES //
////////////////////////////////////////////////////////

// Configuration Variables
struct ConfigVars {
    std::string key;
    std::string value;
};

ConfigVars* InitVar[1024];
int numchar = 0;

void PrintParsingWarning(std::string key)
{
    printf("Warning: No keyword \"%s\" found in config file \n", key.c_str());
}

// Parsing configuration file
void ParserConfigFile(std::string fileName)
{
    std::string optionValue;
    std::ifstream infile;
    infile.open(fileName.c_str());
    numchar = 0;

    // Check if file exist?
    if (infile.is_open() != true)
    {
        printf(" \t Error: %s file does not exist. Quit!!!\n\n", fileName.c_str());
        exit(0);
    }

    std::string key;
    while (!infile.eof())   // While loop to get all the lines.
    {
        getline(infile, optionValue);    // Saves the line in STRING.

        // If the option is a comment, continue
        if (optionValue.substr(0, 1) == "#")
            continue;

        key = ParseName(optionValue);
        if (key.length() > 0)
        {
            InitVar[numchar] = new ConfigVars;
            InitVar[numchar]->key = key;
            InitVar[numchar]->value = ParseValue(optionValue);
            numchar++;
        }
    }

    numchar--;
    infile.close();

}


// Convert Option value to strings
std::string GetOptionToString(std::string key) {
    int flag = 0;

    // Check to see if anything got parsed?
    if (numchar == 0)
        return "";

    for (int i = 0; i <= numchar; i++)
    {
        if (key == InitVar[i]->key)
            return InitVar[i]->value;
    }

    if (flag == 0) 
    {
        PrintParsingWarning(key);
    }    
    return "";
}

// Convert Option value to characters
const char *GetOptionToChar(std::string key)
{
    int flag = 0;

    // Check to see if anything got parsed?
    if (numchar == 0)     
        return "";

    for (int i = 0; i <= numchar; i++)
    {
        if (key == InitVar[i]->key)
            return InitVar[i]->value.c_str();
    }
    
    if (flag == 0) 
    {
        PrintParsingWarning(key);
    }
    return "";
}


// Convert Option value to integer values
int GetOptionToInt(std::string key)
{
    int flag = 0;

    // Check to see if anything got parsed?
    if (numchar == 0)
        return 0;

    for (int i = 0; i <= numchar; i++)
    {
        if (key == InitVar[i]->key)
            return atoi(InitVar[i]->value.c_str());
    }

    if (flag == 0) 
    {
        PrintParsingWarning(key);
    }
    return 0;
}


// Convert Option value to double values
double GetOptionToDouble(std::string key)
{
    int flag = 0;

    // Check to see if anything got parsed?
    if (numchar == 0)
        return 0;

    for (int i = 0; i <= numchar; i++)
    {
        if (key == InitVar[i]->key)
            return atof(InitVar[i]->value.c_str());
    }

    if (flag == 0) 
    {
        PrintParsingWarning(key);
    }

    return 0;
}


// Parse the name
std::string ParseName(std::string value)
{
    size_t str;
    str = value.find('=');

    if (str > 100)
        return "";

    std::string KeyName = value.substr(0, (str-1));
    KeyName = Trim(KeyName);

    return KeyName;
}


// Parse the value
std::string ParseValue(std::string value)
{
    size_t str;
    str = value.find('=');

    if (str > 100)
        return "";

    std::string KeyValue = value.substr((str+1));
    KeyValue = Trim(KeyValue);

    return KeyValue;
}


// trim string
std::string Trim(std::string str)
{
    return LTrim(RTrim(str));
}


// trim from start
std::string LTrim(std::string str)
{
    str.erase(str.begin(), find_if(str.begin(), str.end(), not1(std::ptr_fun<int, int>(isspace))));
    return str;
}


// trim from end
std::string RTrim(std::string str)
{
    str.erase(find_if(str.rbegin(), str.rend(), not1(std::ptr_fun<int, int>(isspace))).base(), str.end());
    return str;
}