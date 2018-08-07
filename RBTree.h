// RBTree.h
#ifndef RBTREE_H_
#define RBTREE_H_

namespace RBTREE
{
	enum COLOR
	{
		RED, BLACK
	};

	template <class T>
	class RBTree
	{
		class RBTree_Node
		{
		private:
			COLOR color;
			T * data;
			RBTree_Node * leftChild;
			RBTree_Node * rightChild;
			RBTree_Node * parent;
			void setParent(RBTree_Node * node)				{ parent = node; }
		public:
			RBTree_Node(const T & dat) : data(new T(dat))	{ leftChild = nullptr; rightChild = nullptr; parent = nullptr; }
			void setLeftChild(RBTree_Node * node)			{ leftChild = node; if(node) node->setParent(this); }
			void setRightChild(RBTree_Node * node)			{ rightChild = node; if (node) node->setParent(this); }
			void clearParent()								{ parent = nullptr; }
			RBTree_Node * getLeftChild()					{ return leftChild; }
			RBTree_Node * getRightChild()					{ return rightChild; }
			RBTree_Node * getParent()						{ return parent; }
			T * getData()									{ return data; }
			COLOR & getColor()								{ return color; }
			RBTree_Node * getUncle();
			RBTree_Node * getGrandParent();
		};

	private:
		RBTree_Node * root;
		bool(*compare)(const T &, const T &);
		bool(*equals)(const T &, const T &);
		void insertHelper(RBTree_Node * root, RBTree_Node * parent, const T & data);
		void insertCase1(RBTree_Node * root);
		void insertCase2(RBTree_Node * root);
		void insertCase3(RBTree_Node * root);
		void insertCase4(RBTree_Node * root);
		void insertCase5(RBTree_Node * root);
		void rebalance(RBTree_Node * root);
		RBTree_Node * searchHelper(RBTree_Node * root, const T & key);
		void inOrderTraversalHelper(RBTree_Node * root, void (*action)( T &));
	public:
		RBTree(const T & rt,
			   bool (*comp)(const T &, const T &) = 
			   [](const T & a, const T & b)
				{
					return (a < b);
				},
				bool (*eq)(const T &, const T &) =
				[](const T & a, const T & b)
				{
					return (a == b);
				})
		: compare(comp), equals(eq)
		{
			insertHelper(nullptr, nullptr, rt);
		}

		// interface presented to the user which 
		// hides the ugly recursive underbelly
		void insert(const T & data);
		T * search(const T & key);
		void inOrderTraversal(void (*action)( T & data));
	};	

	template <class T> typename
	RBTree<T>::RBTree_Node * RBTree<T>::RBTree_Node::getGrandParent()
	{
		if (parent != nullptr)
			return parent->getParent();

		return nullptr;
	}

	template <class T> typename
	RBTree<T>::RBTree_Node * RBTree<T>::RBTree_Node::getUncle()
	{
		RBTree_Node * g = getGrandParent();
		if (g == nullptr)
			return nullptr; // no grandpa means no uncle
		if (parent == g->getLeftChild())
			return g->getRightChild();

		return g->getLeftChild();
	}

	template <class T> typename
	void RBTree<T>::inOrderTraversalHelper(RBTree_Node * root, void(*action)( T &))
	{
		if (!root)
			return;
		inOrderTraversalHelper(root->getLeftChild(), action);
		action(*(root->getData()));
		inOrderTraversalHelper(root->getRightChild(), action);
	}

	template <class T> typename
	void RBTree<T>::inOrderTraversal(void(*action)( T &))
	{
		inOrderTraversalHelper(root, action);
	}

	template <class T> typename
	RBTree<T>::RBTree_Node * RBTree<T>::searchHelper(RBTree_Node * root, const T & key)
	{
		if (!root || equals(*(root->getData()), key))
			return root;
		else if (compare(key, *(root->getData())))
			return searchHelper(root->getLeftChild(), key);
		else
			return searchHelper(root->getRightChild(), key);
	}

	template <class T> typename
	T * RBTree<T>::search(const T & data)
	{
		RBTree_Node * ret = searchHelper(root, data);
		if (ret)
			return (ret->getData());
		return nullptr;
	}

	template <class T> typename
	void RBTree<T>::insertCase5(RBTree_Node * newNode)
	{
		RBTree_Node * g = newNode->getGrandParent();
		
		// is new root of the tree?
		RBTree_Node * temp = root;
		if (g == root)
		{
			temp = newNode->getParent();
		}

		newNode->getParent()->getColor() = BLACK;
		g->getColor() = RED;
		if (newNode == newNode->getParent()->getLeftChild())
		{
			// rotate right on grandparent node
			RBTree_Node * p = newNode->getParent();
			RBTree_Node * x = p->getRightChild();
			RBTree_Node * ggp = g->getParent();
			if (ggp)
			{
				if (ggp->getRightChild() == g)
				{
					ggp->setRightChild(p);
				}
				else
				{
					ggp->setLeftChild(p);
				}
			}
			else
			{
				p->clearParent();
			}
			p->setRightChild(g);
			g->setLeftChild(x);
		}
		else
		{
			// rotate left on grandparent node
			RBTree_Node * p = newNode->getParent();
			RBTree_Node * x = p->getLeftChild();
			RBTree_Node * ggp = g->getParent();
			if (ggp)
			{
				if (ggp->getRightChild() == g)
				{
					ggp->setRightChild(p);
				}
				else
				{
					ggp->setLeftChild(p);
				}
			}
			else
			{
				p->clearParent();
			}
			p->setLeftChild(g);
			g->setRightChild(x);
		}

		// set new root of the tree
		root = temp;
	}

	template <class T> typename
	void RBTree<T>::insertCase4(RBTree_Node * newNode)
	{
		RBTree_Node * g = newNode->getGrandParent();
		if (newNode == newNode->getParent()->getRightChild()
			&& newNode->getParent() == g->getLeftChild())
		{
			// rotate left on parent node
			RBTree_Node * p = newNode->getParent();
			RBTree_Node * left = newNode->getLeftChild();
			g->setLeftChild(newNode);
			newNode->setLeftChild(p);
			p->setRightChild(left);
			newNode = p;
		}
		else if (newNode == newNode->getParent()->getLeftChild()
				 && newNode->getParent() == g->getRightChild())
		{
			// rotate right on parent node
			RBTree_Node * p = newNode->getParent();
			RBTree_Node * right = newNode->getRightChild();
			g->setRightChild(newNode);
			newNode->setRightChild(p);
			p->setLeftChild(right);
			newNode = p;
		}
		insertCase5(newNode);
	}

	template <class T> typename
	void RBTree<T>::insertCase3(RBTree_Node * newNode)
	{
		RBTree_Node * u = newNode->getUncle(),
					* g = newNode->getGrandParent();
		if (u && u->getColor() == RED)
		{
			newNode->getParent()->getColor() = BLACK;
			u->getColor() = BLACK;
			g->getColor() = RED;
			insertCase1(g);
		}
		else
			insertCase4(newNode);
	}

	template <class T> typename
	void RBTree<T>::insertCase2(RBTree_Node * newNode)
	{
		if (newNode->getParent()->getColor() == BLACK)
			return;
		else
			insertCase3(newNode);
	}

	template <class T> typename
	void RBTree<T>::insertCase1(RBTree_Node * newNode)
	{
		if (!(newNode->getParent()))
		{
			root = newNode;
			newNode->getColor() = BLACK;
		}
		else
			insertCase2(newNode);
	}

	template <class T> typename
	void RBTree<T>::rebalance(RBTree_Node * newNode)
	{
		insertCase1(newNode);
	}

	template <class T> typename
	void RBTree<T>::insertHelper(RBTree_Node * root, RBTree_Node * parent, const T & data)
	{
		if (!root)
		{
			root = new RBTree_Node(data);
			if (parent)
				if (compare(data, *(parent->getData())))
					parent->setLeftChild(root);
				else
					parent->setRightChild(root);
			root->getColor() = RED;

			/*	we do RB Tree rebalancing here to ensure
				the tree is correctly balanced
			*/
			rebalance(root);
			
		}
		else if (compare(data, *(root->getData())))
			insertHelper(root->getLeftChild(), root, data);
		else if (compare(*(root->getData()), data))
			insertHelper(root->getRightChild(), root, data);

			// Otherwise ERROR - Reinsertion of the same key
	}

	template <class T> typename
	void RBTree<T>::insert(const T & data)
	{
		if (root)
			insertHelper(root, root->getParent(), data);
		else
			insertHelper(nullptr, nullptr, data);
	}

}


#endif
