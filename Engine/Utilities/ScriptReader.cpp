#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <boost/algorithm/string.hpp>

#include "ScriptReader.h"

void ScriptReader::parseTokens()
{
	using namespace std;

	head_block_ = new ScriptBlock;
	head_block_->name_ = "head_block";

	parseTokens(tokens_, head_block_);
}

void ScriptReader::printAll()
{
	head_block_->print();
}

void ScriptReader::parseTokens(std::list<std::string>& block_tokens, ScriptBlock* block)
{
	using namespace std;

// 	std::cout << "start parsing block " << block->name_ << std::endl;
// 	for (list<string>::iterator itr = block_tokens.begin(); itr != block_tokens.end(); ++itr)
// 	{
// 		std::cout << *itr << " ";
// 	}
// 	std::cout << std::endl;

	// start parsing block tokens
	for (list<string>::iterator itr = block_tokens.begin(); itr != block_tokens.end(); ++itr)
	{
		if (isKeyword(*itr) == true)
		{
			const std::string &word = *itr;

			ScriptParameter par;

			int words_to_read = 1;

			// read type
			if (word == "int") par.type_ = ScriptParameter::Type::INT;
			else if (word == "int2")
			{
				par.type_ = ScriptParameter::Type::INT2;
				words_to_read = 2;
			}
			else if (word == "int3")
			{
				par.type_ = ScriptParameter::Type::INT3;
				words_to_read = 3;
			}
			else if (word == "bool") par.type_ = ScriptParameter::Type::BOOL;
			else if (word == "float") par.type_ = ScriptParameter::Type::FLOAT;
			else if (word == "double") par.type_ = ScriptParameter::Type::DOUBLE;
			else if (word == "string") par.type_ = ScriptParameter::Type::STRING;
			else if (word == "float2")
			{
				par.type_ = ScriptParameter::Type::FLOAT2;
				words_to_read = 2;
			}
			else if (word == "float3")
			{
				par.type_ = ScriptParameter::Type::FLOAT3;
				words_to_read = 3;
			}
			else if (word == "float4")
			{
				par.type_ = ScriptParameter::Type::FLOAT4;
				words_to_read = 4;
			}
			else
			{
				std::cout << "Wrong script parameter type definition " << word << std::endl;
				exit(1);
			}

			itr++; if (itr == block_tokens.end()){ std::cout << "Unexpected end of tokens " << std::endl; exit(1); }

			// read name
			par.name_ = *itr;

			// read values
			for (int w = 0; w < words_to_read; w++) {

				itr++; if (itr == block_tokens.end()){ std::cout << "Unexpected end of tokens " << std::endl; exit(1); }

				std::string read_block = *itr;

				par.value_.append(read_block);

				if (w != words_to_read - 1) par.value_.append(" ");	// do not add blank at the end of value_ string.
			}

			block->par_list_.par_list_.push_back(par);
		}
		else// if token is not one of the keywords, it is a block name
		{
			ScriptBlock *child_block = new ScriptBlock;

			block->children_.push_back(child_block);

			child_block->name_ = *itr;	// block name

			itr++;

			if (itr == block_tokens.end()){ std::cout << "Unexpected script end" << std::endl; exit(1); }

			list<string> sub_tokens;

			// find tokens in the block
			if (boost::istarts_with(*itr, "{") == true)	// if block starts
			{
				itr++;	// skip {

				if (itr == block_tokens.end()){ std::cout << "Unexpected script end" << std::endl; exit(1); }

				int bracket_count = 0;

				while (true)
				{
					if (boost::istarts_with(*itr, "{") == true)
					{
						bracket_count++;
					}

					if (boost::istarts_with(*itr, "}") == true)
					{
						if (bracket_count == 0) break; 
						
						bracket_count--;

						if (bracket_count < 0){ std::cout << "Wrong bracket count " << std::endl; exit(1); }
					}

					sub_tokens.push_back(*itr);

					itr++;

					if (itr == block_tokens.end())
					{
						assert(false);
						std::cout << "Matching } was not found" << std::endl; exit(1);
					}

				}

				parseTokens(sub_tokens, child_block);
			}
			else
			{
				std::cout << "Can't find block start {" << std::endl;
				exit(1);
			}
		}
	}
}

bool ScriptReader::isKeyword(const std::string& input)
{
	using namespace std;

	for (list<string>::iterator itr = keywords_.begin(); itr != keywords_.end(); ++itr)
	{
//		if (boost::istarts_with(input, *itr) == true) return true;

		if ((*itr).compare(input) == 0) return true;
	}

	return false;
}

void ScriptReader::findAllTokens(std::ifstream& file)
{
	using namespace std;

	std::vector<std::string> block;
	string input_line;
	bool comment_block = false;
	while (getline(file, input_line))		// line loop
	{
		using namespace std;

		std::istringstream iss(input_line, std::istringstream::in);

		string word;

		while (iss >> word)	// words loop in the line
		{
			if (boost::istarts_with(word, "//") || boost::istarts_with(word, "#")) break;		// skip comment lines

			if (comment_block == false && boost::istarts_with(word, "/*"))		// start comment block
			{
				comment_block = true;
				continue;
			}
			else if (comment_block == true && boost::istarts_with(word, "*/"))		// end comment block
			{
				comment_block = false;
				boost::replace_first(word, "*/", "");
			}
			else if (comment_block == true)	// ignore lines in comment block
			{
				continue;
			}

			if(!word.empty()) tokens_.push_back(word);
		}
	}
}

void ScriptReader::readFile(const char* filename)
{
	using namespace std;

	std::cout << "Start reading " << filename << std::endl;

	ifstream file(filename);

	// check if file is correctly opened
	if (file.is_open() == false){ std::cout << filename << " does not exist. Program terminated." << std::endl; exit(-1); }

	findAllTokens(file);

	file.close();

	// printAllTokens
// 	for (list<string>::iterator itr = tokens_.begin(); itr != tokens_.end(); ++itr)
// 		std::cout << *itr << " ";
// 	std::cout << std::endl;

	parseTokens();

//	printAll();
}