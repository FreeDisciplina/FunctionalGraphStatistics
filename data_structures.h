#ifndef DATA_STRUCTURES_H__
#define DATA_STRUCTURES_H__

#include "globalvars.h"

typedef struct fg_node
{
	data_t image;
	u64 preimageN;
	u64 depth;
} info_node;

typedef struct BST
{
	data_t val;
	struct BST * left;
	struct BST * right;
	BST() :left(NULL), right(NULL) {};
} BSTree;

typedef struct BST2
{
	data_t val;
	data_t num;
	struct BST2 * left;
	struct BST2 * right;
	BST2() :left(NULL), right(NULL) {};
} BSTree2;

typedef struct T
{
	data_t val;
	struct T * next;
	T() :next(NULL) {};
} T_node;


typedef struct list
{
	data_t val;
	struct list *next;
	list() :next(NULL) {};
} preimageList;

T_node * addToT(T_node * node, data_t val)
{
	T_node * tmp;
	tmp = (T_node *)malloc(sizeof(T_node));
	if (tmp == NULL)
	{
		cout << "Memory allocation error. T_node * tmp is NULL!" << endl;
#if defined(_MSC_VER)
		system("Pause");
#endif
		exit(1);
	}

	tmp->val = val;
	tmp->next = node;
	return tmp;
}

void clearT(T_node *node)
{
	while (node != NULL)
	{
		T_node * tmp = node;
		node = node->next;
		tmp->next = NULL;
		free(tmp);
	}
}


BSTree * getNewNode(data_t val)
{
	BSTree * tmp = (BSTree *)malloc(sizeof(BSTree));
	if (tmp == NULL)
	{
		cout << "Memory allocation error. BSTree * tmp is NULL!" << endl;
#if defined(_MSC_VER)
		system("Pause");
#endif
		exit(1);
	}

	tmp->val = val;
	tmp->left = NULL;
	tmp->right = NULL;
	return tmp;
}

BSTree * insertBST(bool * isNotIn, BSTree * root, data_t val)
{
	if (root == NULL)
	{
		root = getNewNode(val);
		*isNotIn = true;
	}
	else if (val == root->val)
	{
		*isNotIn = false;
		return root;
	}
	else if (val < root->val)
	{
		root->left = insertBST(isNotIn, root->left, val);
	}
	else
	{
		root->right = insertBST(isNotIn, root->right, val);
	}
	return root;
}

BSTree * extractBST(BSTree * root, data_t val)
{
	if (root == NULL)
	{
		return root;
	}
	else if (val < root->val)
	{
		root->left = extractBST(root->left, val);
	}
	else if (val > root->val)
	{
		root->right = extractBST(root->right, val);
	}
	else
	{
		if ((root->left == NULL) && (root->right == NULL))
		{
			BSTree * tmp = root;
			root = NULL;
			free(tmp);
		}
		else if (root->left == NULL)
		{
			BSTree * tmp = root;
			root = root->right;
			tmp->right = NULL;
			free(tmp);
		}
		else if (root->right == NULL)
		{
			BSTree * tmp = root;
			root = root->left;
			tmp->left = NULL;
			free(tmp);
		}
		else
		{
			BSTree * tmp = root->right;
			while (tmp->left != NULL) tmp = tmp->left;
			root->val = tmp->val;
			root->right = extractBST(root->right, tmp->val);
		}
	}
	return root;
}

void clearBST(BSTree * root)
{
	if (root == NULL) return;
	clearBST(root->left);
	clearBST(root->right);
	free(root);
}


BSTree2 * getNewNodeBSTree2(data_t val)
{
	BSTree2 * tmp = (BSTree2 *)malloc(sizeof(BSTree2));
	if (tmp == NULL)
	{
		cout << "Memory allocation error. BSTree2 * tmp is NULL!" << endl;
#if defined(_MSC_VER)
		system("Pause");
#endif
		exit(1);
	}

	tmp->val = val;
	tmp->num = 1;
	tmp->left = NULL;
	tmp->right = NULL;
	return tmp;
}

BSTree2 * insertBST2(bool * isNotIn, BSTree2 * root, data_t val)
{
	if (root == NULL)
	{
		root = getNewNodeBSTree2(val);
		*isNotIn = true;
	}
	else if (val == root->val)
	{
		*isNotIn = false;
		root->num++;
		return root;
	}
	else if (val < root->val)
	{
		root->left = insertBST2(isNotIn, root->left, val);
	}
	else
	{
		root->right = insertBST2(isNotIn, root->right, val);
	}
	return root;
}

BSTree2 * extractBST2(BSTree2 * root, data_t val)
{
	if (root == NULL)
	{
		return root;
	}
	else if (val < root->val)
	{
		root->left = extractBST2(root->left, val);
	}
	else if (val > root->val)
	{
		root->right = extractBST2(root->right, val);
	}
	else
	{
		if ((root->left == NULL) && (root->right == NULL))
		{
			BSTree2 * tmp = root;
			root = NULL;
			free(tmp);
		}
		else if (root->left == NULL)
		{
			BSTree2 * tmp = root;
			root = root->right;
			tmp->right = NULL;
			free(tmp);
		}
		else if (root->right == NULL)
		{
			BSTree2 * tmp = root;
			root = root->left;
			tmp->left = NULL;
			free(tmp);
		}
		else
		{
			BSTree2 * tmp = root->right;
			while (tmp->left != NULL) tmp = tmp->left;
			root->val = tmp->val;
			root->num = tmp->num;
			root->right = extractBST2(root->right, tmp->val);
		}
	}
	return root;
}

void clearBST2(BSTree2 * root)
{
	if (root == NULL) return;
	clearBST2(root->left);
	clearBST2(root->right);
	free(root);
}

// x is an array of elements that are of type u64.
#define getBit(x, i) (((x)[(i) >> 6] >> ((i) & 0x3fULL)) & 0x1ULL)
#define setBit(x, i) x[(i) >> 6] = x[(i) >> 6] & ((0x1ULL << ((i) & 0x3fULL)) ^ 0xffffffffffffffffULL);

u64 getIndicated(data_t * nodes, u64 * indicator)
{
	u64 number = 0ULL;
	u64 base = 0ULL;
	for (u64 i = 0; i < (N >> 6ULL); i++)
	{
		u64 ind = indicator[i];
		u64 indv;
		u64 wt;
		__m256i offtmp0;
		__m256i offtmp1;
		__m256i offtmp01;
		data_t offtmp[16];

		indv = ind & 0x000000000000ffffULL;
		wt = (u64)_popcnt64(indv) & 0x00000000ffffffffULL;
		if (wt != 0ULL)
		{
			offtmp0 = _mm256_loadu_si256((__m256i *)(W16v[indv]));
			offtmp1 = _mm256_loadu_si256((__m256i *)(&(W16v[indv][8])));
			offtmp01 = _mm256_set1_epi32((u32)base);
			offtmp0 = _mm256_add_epi32(offtmp0, offtmp01);
			offtmp1 = _mm256_add_epi32(offtmp1, offtmp01);
			_mm256_storeu_si256((__m256i *)(&(offtmp[0])), offtmp0);
			_mm256_storeu_si256((__m256i *)(&(offtmp[8])), offtmp1);
			std::memcpy((void *)(&(nodes[number])), offtmp, wt * sizeof(data_t));
			number += wt;
		}

		base += 16;
		indv = (ind >> 16ULL) & 0x000000000000ffffULL;
		wt = (u64)_popcnt64(indv) & 0x00000000ffffffffULL;
		if (wt != 0ULL)
		{
			offtmp0 = _mm256_loadu_si256((__m256i *)(W16v[indv]));
			offtmp1 = _mm256_loadu_si256((__m256i *)(&(W16v[indv][8])));
			offtmp01 = _mm256_set1_epi32((u32)base);
			offtmp0 = _mm256_add_epi32(offtmp0, offtmp01);
			offtmp1 = _mm256_add_epi32(offtmp1, offtmp01);
			_mm256_storeu_si256((__m256i *)(&(offtmp[0])), offtmp0);
			_mm256_storeu_si256((__m256i *)(&(offtmp[8])), offtmp1);
			std::memcpy((void *)(&(nodes[number])), offtmp, wt * sizeof(data_t));
			number += wt;
		}

		base += 16;
		indv = (ind >> 32ULL) & 0x000000000000ffffULL;
		wt = (u64)_popcnt64(indv) & 0x00000000ffffffffULL;
		if (wt != 0ULL)
		{
			offtmp0 = _mm256_loadu_si256((__m256i *)(W16v[indv]));
			offtmp1 = _mm256_loadu_si256((__m256i *)(&(W16v[indv][8])));
			offtmp01 = _mm256_set1_epi32((u32)base);
			offtmp0 = _mm256_add_epi32(offtmp0, offtmp01);
			offtmp1 = _mm256_add_epi32(offtmp1, offtmp01);
			_mm256_storeu_si256((__m256i *)(&(offtmp[0])), offtmp0);
			_mm256_storeu_si256((__m256i *)(&(offtmp[8])), offtmp1);
			std::memcpy((void *)(&(nodes[number])), offtmp, wt * sizeof(data_t));
			number += wt;
		}

		base += 16;
		indv = (ind >> 48ULL) & 0x000000000000ffffULL;
		wt = (u64)_popcnt64(indv) & 0x00000000ffffffffULL;
		if (wt != 0ULL)
		{
			offtmp0 = _mm256_loadu_si256((__m256i *)(W16v[indv]));
			offtmp1 = _mm256_loadu_si256((__m256i *)(&(W16v[indv][8])));
			offtmp01 = _mm256_set1_epi32((u32)base);
			offtmp0 = _mm256_add_epi32(offtmp0, offtmp01);
			offtmp1 = _mm256_add_epi32(offtmp1, offtmp01);
			_mm256_storeu_si256((__m256i *)(&(offtmp[0])), offtmp0);
			_mm256_storeu_si256((__m256i *)(&(offtmp[8])), offtmp1);
			std::memcpy((void *)(&(nodes[number])), offtmp, wt * sizeof(data_t));
			number += wt;
		}

		base += 16;
	}
	return number;
}

u64 getIndicatedNumber(u64 * indicator)
{
	u64 number = 0ULL;
	for (u64 i = 0; i < (N >> 6); i++)
	{
		number += (u64)_popcnt64(indicator[i]) & 0x00000000ffffffffULL;
	}
	return number;
}

BSTree * getIndicated_cir(u64 * getNumber, BSTree * root, u64 * indicator)
{
	u64 number = 0ULL;
	bool isNotIn;
	u64 base = 0x0ULL;
	for (u64 i = 0; i < (N >> 6); i++)
	{
		u64 ind = indicator[i];
		u64 indv;
		u64 wt;
		__m256i offtmp0;
		__m256i offtmp1;
		__m256i offtmp01;
		data_t offtmp[16];

		indv = ind & 0x000000000000ffffULL;
		wt = (u64)_popcnt64(indv) & 0x00000000ffffffffULL;
		if (wt != 0)
		{
			offtmp0 = _mm256_loadu_si256((__m256i *)(W16v[indv]));
			offtmp1 = _mm256_loadu_si256((__m256i *)(&(W16v[indv][8])));
			offtmp01 = _mm256_set1_epi32((u32)base);
			offtmp0 = _mm256_add_epi32(offtmp0, offtmp01);
			offtmp1 = _mm256_add_epi32(offtmp1, offtmp01);
			_mm256_storeu_si256((__m256i *)&(offtmp[0]), offtmp0);
			_mm256_storeu_si256((__m256i *)&(offtmp[8]), offtmp1);
			for (u64 j = 0; j < wt; j++)
			{
				root = insertBST(&isNotIn, root, offtmp[j]);
			}
			number += wt;
		}

		base += 16;
		indv = (ind >> 16) & 0x000000000000ffffULL;
		wt = (u64)_popcnt64(indv) & 0x00000000ffffffffULL;
		if (wt != 0)
		{
			offtmp0 = _mm256_loadu_si256((__m256i *)(W16v[indv]));
			offtmp1 = _mm256_loadu_si256((__m256i *)(&(W16v[indv][8])));
			offtmp01 = _mm256_set1_epi32((u32)base);
			offtmp0 = _mm256_add_epi32(offtmp0, offtmp01);
			offtmp1 = _mm256_add_epi32(offtmp1, offtmp01);
			_mm256_storeu_si256((__m256i *)&(offtmp[0]), offtmp0);
			_mm256_storeu_si256((__m256i *)&(offtmp[8]), offtmp1);
			for (u64 j = 0; j < wt; j++)
			{
				root = insertBST(&isNotIn, root, offtmp[j]);
			}
			number += wt;
		}

		base += 16;
		indv = (ind >> 32) & 0x000000000000ffffULL;
		wt = (u64)_popcnt64(indv) & 0x00000000ffffffffULL;
		if (wt != 0)
		{
			offtmp0 = _mm256_loadu_si256((__m256i *)(W16v[indv]));
			offtmp1 = _mm256_loadu_si256((__m256i *)(&(W16v[indv][8])));
			offtmp01 = _mm256_set1_epi32((u32)base);
			offtmp0 = _mm256_add_epi32(offtmp0, offtmp01);
			offtmp1 = _mm256_add_epi32(offtmp1, offtmp01);
			_mm256_storeu_si256((__m256i *)&(offtmp[0]), offtmp0);
			_mm256_storeu_si256((__m256i *)&(offtmp[8]), offtmp1);
			for (u64 j = 0; j < wt; j++)
			{
				root = insertBST(&isNotIn, root, offtmp[j]);
			}
			number += wt;
		}

		base += 16;
		indv = (ind >> 48) & 0x000000000000ffffULL;
		wt = (u64)_popcnt64(indv) & 0x00000000ffffffffULL;
		if (wt != 0)
		{
			offtmp0 = _mm256_loadu_si256((__m256i *)(W16v[indv]));
			offtmp1 = _mm256_loadu_si256((__m256i *)(&(W16v[indv][8])));
			offtmp01 = _mm256_set1_epi32((u32)base);
			offtmp0 = _mm256_add_epi32(offtmp0, offtmp01);
			offtmp1 = _mm256_add_epi32(offtmp1, offtmp01);
			_mm256_storeu_si256((__m256i *)&(offtmp[0]), offtmp0);
			_mm256_storeu_si256((__m256i *)&(offtmp[8]), offtmp1);
			for (u64 j = 0; j < wt; j++)
			{
				root = insertBST(&isNotIn, root, offtmp[j]);
			}
			number += wt;
		}

		base += 16;
	}
	*getNumber = number;
	return root;
}

#endif
