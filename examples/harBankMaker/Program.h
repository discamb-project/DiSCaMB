#pragma once

#include "json.hpp"

#include <string>
#include <map>

class Program 
{
public: 
    virtual ~Program()=0;
    
    virtual void set() = 0;
    
    virtual void run() = 0;

    static Program* create(const std::string& name);
    static void readSettings(nlohmann::json &settings);
};

