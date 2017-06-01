#include "ScriptBlock.h"
#include <sstream>

ScriptBlock::ScriptBlock()
	: parent_(nullptr), name_(std::string("noname"))
{
}

ScriptBlock::~ScriptBlock()
{
	for (std::list<ScriptBlock*>::iterator itr = children_.begin(); itr != children_.end(); ++itr)
	{
		delete (*itr);
	}
}
