#include "ScriptParameterList.h"
#include "exprtk/exprtk.hpp"
#include <sstream>

ScriptParameterList::ScriptParameterList()
	: use_cout_(false)
{}

ScriptParameterList::ScriptParameterList(const char* filename, const bool _use_cout)
	: use_cout_(_use_cout)
{
	initializeFromFile(filename);
}

void ScriptParameterList::initialize(const char* filename, const bool _use_cout)
{
	use_cout_ = _use_cout;

	initializeFromFile(filename);
}

void ScriptParameterList::initializeFromFile(const char *filename)
{
	using namespace std;

	if (use_cout_) std::cout << "Start reading parameter file " << filename << std::endl;

	vector<string> block;

	initializeBlockFromFile(filename, block);

	parseBlock(block);

	if (use_cout_ == true) printAll();
}

void ScriptParameterList::initializeBlockFromFile(const char *filename, std::vector<std::string>& block)
{
	using namespace std;

	ifstream file(filename);

	// check if file is correctly opened
	if (file.is_open() == false){ std::cout << filename << " does not exist. Program terminated." << std::endl; exit(-1); }

	while (true)
	{
		string line;

		if (!getline(file, line)) break;

		block.push_back(line);
	}

	file.close();
}

void ScriptParameterList::initializeFromArgumentsAndFile(const int argc, char* argv[], const std::string default_script_filename)
{
	ScriptParameterList arg_reader;
	arg_reader.initializeFromArguments(argc, argv);

	const std::string script_filename = arg_reader.isValid(std::string("script")) ? arg_reader.getStringValue(std::string("script")) : default_script_filename;

	initializeFromFile(script_filename.c_str());
	combineWith(arg_reader);
}

void ScriptParameterList::printAll()
{
	std::cout << "Print parameter list" << std::endl;

	for (std::list<ScriptParameter>::iterator itr = par_list_.begin(); itr != par_list_.end(); ++itr)
		std::cout << itr->name_ << " " << itr->value_ << std::endl;
}

void ScriptParameterList::parseBlock(const std::vector<std::string>& block)
{
	for (int i = 0; i < (int)block.size(); i++) parseLine(block[i]);
}

void ScriptParameterList::parseTokens(std::list<std::string>& tokens)
{
	for (std::list<std::string>::iterator itr = tokens.begin(); itr != tokens.end(); ++itr)
	{
		const std::string &word = *itr;

		ScriptParameter par;

		int words_to_read = 1;

		// read type
		if (word == "int") par.type_ = ScriptParameter::Type::INT;
		else if (word == "bool") par.type_ = ScriptParameter::Type::BOOL;
		else if (word == "float") par.type_ = ScriptParameter::Type::FLOAT;
		else if (word == "double") par.type_ = ScriptParameter::Type::DOUBLE;
		else if (word == "string") par.type_ = ScriptParameter::Type::STRING;
		else if (word == "float4"){
			par.type_ = ScriptParameter::Type::FLOAT4;
			words_to_read = 4;
		}
		else if (word == "float2"){
			par.type_ = ScriptParameter::Type::FLOAT2;
			words_to_read = 2;
		}
		else {
			std::cout << "Wrong script parameter type definition " << word << std::endl;
			exit(1);
		}

		itr++;if (itr == tokens.end()){ std::cout << "Unexpected end of tokens " << std::endl; exit(1); }
		
		// read name
		par.name_ = *itr;

		// read values
		for (int w = 0; w < words_to_read; w++) {

			itr++; if (itr == tokens.end()){ std::cout << "Unexpected end of tokens " << std::endl; exit(1); }

			std::string read_block = *itr;

			par.value_.append(read_block);

			if (w != words_to_read - 1) par.value_.append(" ");	// do not add blank at the end of value_ string.
		}

		par_list_.push_back(par);
	}
}

void ScriptParameterList::parseLine(const std::string& line)
{
	using namespace std;

	istringstream iss(line, istringstream::in);

	string word;

	// TODO: ignore comment block
// 		if (boost::istarts_with(word, "/*")) // skip comments block
// 			while (iss >> word){
// 				if (boost::iends_with(word, "*/"))break;
// 			};

	while (iss >> word)
	{
		if (boost::istarts_with(word, "//") || boost::istarts_with(word, "#")) break;		// skip comments

		ScriptParameter par;

		int words_to_read = 1;

		// read type
		if (word == "int") par.type_ = ScriptParameter::Type::INT;
		else if (word == "bool") par.type_ = ScriptParameter::Type::BOOL;
		else if (word == "float") par.type_ = ScriptParameter::Type::FLOAT;
		else if (word == "double") par.type_ = ScriptParameter::Type::DOUBLE;
		else if (word == "string") par.type_ = ScriptParameter::Type::STRING;
		else if (word == "float4"){
			par.type_ = ScriptParameter::Type::FLOAT4;
			words_to_read = 4;
		}
		else if (word == "float2"){
			par.type_ = ScriptParameter::Type::FLOAT2;
			words_to_read = 2;
		}
		else {
			std::cout << "Wrong script parameter type definition " << word << std::endl;
			exit(1);
		}

		// read name
		iss >> par.name_;

		// read values
		for (int w = 0; w < words_to_read; w++) {
			string read_block;
			iss >> read_block;

			par.value_.append(read_block);

			if (w != words_to_read - 1) par.value_.append(" ");	// do not add blank at the end of value_ string.
		}

		par_list_.push_back(par);
	}
}

void ScriptParameterList::setStringValue(const std::string par_name, const std::string par_value)
{
	for (std::list<ScriptParameter>::iterator itr = par_list_.begin(); itr != par_list_.end(); ++itr)
	{
		if (itr->name_.compare(par_name) == 0) itr->value_ = par_value;
	}
}

std::string ScriptParameterList::getStringValue(const std::string par_name)
{
	for (std::list<ScriptParameter>::iterator itr = par_list_.begin(); itr != par_list_.end(); ++itr)
	{
		if (itr->name_.compare(par_name) == 0) return itr->value_;
	}

	std::cout << par_name << " parameter does not exist" << std::endl;
	exit(1);
}

bool ScriptParameterList::isValid(const std::string par_name)
{
	for (std::list<ScriptParameter>::iterator itr = par_list_.begin(); itr != par_list_.end(); ++itr)
	{
		if (itr->name_.compare(par_name) == 0) return true;
	}

	return false;
}

void ScriptParameterList::combineWith(ScriptParameterList& input_list)
{
	for (std::list<ScriptParameter>::iterator itr = input_list.par_list_.begin(); itr != input_list.par_list_.end(); ++itr)
	{
		if (isValid(itr->name_) == true) setStringValue(itr->name_, itr->value_);		// overwrite input_list parameters
		else par_list_.push_front(*itr);
	}
}

void ScriptParameterList::overwriteFromArtuments(const int argc, char* argv[])
{
}

void ScriptParameterList::initializeFromMessage(const std::string message)
{
	// parse command messages : -type FLUID_SIMULATION -script_filename script.sim -repeated_run true ...

	using namespace std;

	istringstream iss(message, istringstream::in);

	string word;

	string name, value;

	ScriptParameter par;

	while (iss >> word)
	{
		if (boost::istarts_with(word, "-"))
		{
			if (boost::empty(value) == false)
			{
				boost::replace_last(value, " ", "");
				par.type_ = ScriptParameter::Type::STRING;
				par.value_ = value;
				par.name_ = name;
				par_list_.push_back(par);
				value = "";
			}

			name = word;
			boost::replace_first(name, "-", "");
		}
		else
		{
			value += word + " ";
		}
	}

	if (boost::empty(value) == false)
	{
		boost::replace_last(value, " ", "");
		par.type_ = ScriptParameter::Type::STRING;
		par.value_ = value;
		par.name_ = name;
		par_list_.push_back(par);
		value = "";
	}
}

void ScriptParameterList::initializeFromArguments(const int argc, char *argv[])
{
	par_list_.push_back(ScriptParameter(ScriptParameter::Type::STRING, "name_of_execution", argv[0])); // name of execution can be used to find execution directory

//	par_list_.reserve((argc - 1) / 2);		// don't need reserve linked list

	if ((argc - 1) % 2 != 0)	// argc should be a odd number. arguments are : name of execution -flag1 value1 -flag2 value2 ...
	{
		std::cout << "Wrong arguments" << std::endl;
		std::cin.get();
		exit(1);
	}

	// read flag (-name) and value pairs
	int i = 1;
	while (i < argc)
	{
		std::string flag = argv[i++];
		const std::string value = argv[i++];

		if (!boost::istarts_with(flag, "-"))
		{
			std::cout << "Flag " << flag << " does not start with -" << std::endl;
			exit(1);
		}

		boost::replace_first(flag, "-", "");	// remove "-" from flag to obtain parameter name. flag is "-name".

		par_list_.push_back(ScriptParameter(ScriptParameter::Type::UNKNOWN, flag, value));
	}
}

bool ScriptParameterList::getValue(const std::string par_name, const bool type_specifier)
{
	std::string value = getStringValue(par_name);

	if (value.compare("false") == 0 || value.compare("0") == 0) return false;
	else return true;
}

int ScriptParameterList::getValue(const std::string par_name, const int type_specifier)
{
	return (int)getValueFromString(getStringValue(par_name));
}

std::string ScriptParameterList::getValue(const std::string par_name, const std::string type_specifier)
{
	return getStringValue(par_name);
}

float ScriptParameterList::getValue(const std::string par_name, const float type_specifier)
{
	return getValueFromString(getStringValue(par_name));
}


Vector2D<int> ScriptParameterList::getValue(const std::string par_name, const Vector2D<int> type_specifier)
{
	std::vector<std::string> value_splits;

	std::string value = getStringValue(par_name);

	boost::split(value_splits, value, boost::is_any_of(" "), boost::token_compress_on);

	return Vector2D<int>((int)getValueFromString(value_splits[0]), (int)getValueFromString(value_splits[1]));
}

Vector3D<int> ScriptParameterList::getValue(const std::string par_name, const Vector3D<int> type_specifier)
{
	std::vector<std::string> value_splits;

	std::string value = getStringValue(par_name);

	boost::split(value_splits, value, boost::is_any_of(" "), boost::token_compress_on);

	return Vector3D<int>(getValueFromString(value_splits[0]), getValueFromString(value_splits[1]), getValueFromString(value_splits[2]));
}

Vector4D<int> ScriptParameterList::getValue(const std::string par_name, const Vector4D<int> type_specifier)
{
	std::vector<std::string> value_splits;

	std::string value = getStringValue(par_name);

	boost::split(value_splits, value, boost::is_any_of(" "), boost::token_compress_on);

	return Vector4D<int>(getValueFromString(value_splits[0]), getValueFromString(value_splits[1]), getValueFromString(value_splits[2]), getValueFromString(value_splits[3]));
}

Vector2D<float> ScriptParameterList::getValue(const std::string par_name, const Vector2D<float> type_specifier)
{
	std::vector<std::string> value_splits;

	std::string value = getStringValue(par_name);

	boost::split(value_splits, value, boost::is_any_of(" "), boost::token_compress_on);

	return Vector2D<float>(getValueFromString(value_splits[0]), getValueFromString(value_splits[1]));
}

Vector3D<float> ScriptParameterList::getValue(const std::string par_name, const Vector3D<float> type_specifier)
{
	std::vector<std::string> value_splits;

	std::string value = getStringValue(par_name);

	boost::split(value_splits, value, boost::is_any_of(" "), boost::token_compress_on);

	return Vector3D<float>(getValueFromString(value_splits[0]), getValueFromString(value_splits[1]), getValueFromString(value_splits[2]));
}

Vector4D<float> ScriptParameterList::getValue(const std::string par_name, const Vector4D<float> type_specifier)
{
	std::vector<std::string> value_splits;

	std::string value = getStringValue(par_name);

	boost::split(value_splits, value, boost::is_any_of(" "), boost::token_compress_on);

	return Vector4D<float>(getValueFromString(value_splits[0]), getValueFromString(value_splits[1]), getValueFromString(value_splits[2]), getValueFromString(value_splits[3]));
}

float ScriptParameterList::getValueFromString(const std::string string_expression) const
{
	exprtk::expression<float> expression;
	exprtk::parser<float>     parser;

	parser.compile(string_expression, expression);

	return expression.value();
}