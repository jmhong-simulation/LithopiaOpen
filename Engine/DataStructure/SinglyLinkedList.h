// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

/* How to use

for (object_list_.begin(); object_list_.valid();)
{
	Object *temp = object_list_.getItr();
	object_list_.next();
	SAFE_DELETE(temp);
}

object_list_.reset();

*/

template<class TT>
class SinglyLinkedElement
{
public:
	typedef typename SinglyLinkedElement<TT> T_ELEMENT;

	TT value_;
	T_ELEMENT *next_;

	SinglyLinkedElement(const TT& value)
		: value_(value), next_(nullptr)
	{}
};

template<class TT>
class SinglyLinkedList
{
public:
	typedef typename TT TT;
	typedef typename SinglyLinkedElement<TT> T_ELEMENT;

	T_ELEMENT *head_;
	T_ELEMENT *last_;
	T_ELEMENT *itr_;

	SinglyLinkedList()
		: head_(nullptr), last_(nullptr), itr_(nullptr)
	{}

	~SinglyLinkedList()
	{
		reset();
	}

	void reset()
	{
		for (T_ELEMENT *itr = head_; itr != nullptr;)
		{
			T_ELEMENT *temp = itr;

			itr = itr->next_;

			delete temp;
		}

		head_ = nullptr;
		last_ = nullptr;
	}

	void pushFront(const TT value)
	{
		T_ELEMENT *new_element = new T_ELEMENT(value);

		if (head_ == nullptr)
		{
			head_ = new_element;
			last_ = new_element;
		}
		else
		{
			new_element->next_ = head_;
			head_ = new_element;
		}
	}

	void pushBack(const TT value)
	{
		T_ELEMENT *new_element = new T_ELEMENT(value);

		if (last_ == nullptr)
		{
			head_ = new_element;
			last_ = new_element;
		}
		else
		{
			last_->next_ = new_element;
			last_ = new_element;
		}
	}

	// iterator functions
	void begin()
	{
		itr_ = head_;
	}

	bool valid()
	{
		if (itr_ != nullptr) return true;
		else return false;
	}

	void next()
	{
		itr_ = itr_->next_;
	}

	TT& getItr()
	{
		return itr_->value_;
	}

	TT& getItr(const int ix)
	{
		begin();

		for (int i = 0; i < ix; i++)
		{
			assert(valid());

			next();
		}

		return itr_->value_;
	}

	void resetPointers()		// works only when TT is pointer
	{
		for (begin(); valid();)
		{
			TT temp = getItr();
			next();
			SAFE_DELETE(temp);
		}

		reset();
	}
};

