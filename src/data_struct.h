/*
 *	Author: Cheng Long
 *	Email: clong@cse.ust.hk
 */

#ifndef DATA_STRUCT_H
#define DATA_STRUCT_H

//#include "bst.h"

//#define	KEY_TYPE	float
#define	KEY_TYPE	double
#define B_KEY_TYPE	double
//#define KEY2_TYPE	double
#define W_TYPE int

//#define max(a, b) ((a) >= (b) ? (a) : (b))
//#define min(a, b) ((a) <= (b) ? (a) : (b))



//#define	COST_TAG	2

//From costenum.h
//Inverted File related structures.
typedef struct k_node 
{
	KEY_TYPE		key;
	struct k_node*	next;

    int freq;		//unified pruning
	W_TYPE			weight; //IR-tree
	
}	k_node_t;


//The structure for storing a set of keywords.
typedef struct psi
{
    int			key_n;
    k_node_t*	k_head;
    
}	psi_t;


//From irtree.h
//The interval structure of one dimension.
typedef struct 
{
	float min;
	float max;

}	range;

//The structure of the object.
typedef struct obj
{
	int			id;
	range*		MBR;

	//IR-tree augmentation.
	k_node_t*	k_head;
	
	W_TYPE	weight;
	
}	obj_t;

//The node structure for storing an object pointer.
typedef struct obj_node
{
	obj_t*				obj_v;
	struct obj_node*	next;

    psi_t*              psi_v;      ///storing the intersection of query keyword and obj keyword
    
	//for const_k_NN_key only.
	B_KEY_TYPE			dist;
    
    ///for cost 2.
    int                 numOfInt;
    
    //for SKECa+
    W_TYPE          maxInvalidRange;
    bool                type; //true=in, false = out
    B_KEY_TYPE          angle;
    struct obj_set*     range;
	
	struct bst*				IF_S;
	
}	obj_node_t;

//The structure for storing a group of objects.
typedef struct obj_set
{
	int			obj_n;
	obj_node_t*	head;

}	obj_set_t;



///---
typedef struct obj_dist_node
{
	B_KEY_TYPE*				s;
    obj_node_t**             n;
	struct obj_dist_node*	next;

    
}	obj_dist_node_t;

typedef struct obj_dist_set
{
	int			obj_n;
	obj_dist_node_t*	head;
}	obj_dist_set_t;


///-------

//From costenum.h.

//The structure for storing a location.
typedef struct loc
{
	int		dim;
	float*	coord;
}	loc_t;

//Structure for storing a disk.
typedef struct disk
{
	loc_t*		center;
	B_KEY_TYPE	radius;
}	disk_t;


//Structure for storing an ellipse.
typedef struct ellipse
{
	loc_t*		center1;
	loc_t*		center2;
	B_KEY_TYPE	diam;
	
}	ellipse_t;


//Structure for storing a query.
typedef struct query
{
	//loc.
	loc_t*	loc_v;

	//doc.
	//int			key_n;
	//k_node_t*	doc_v;

	psi_t*	psi_v;
}	query_t;

//Triplet structure.
typedef struct tri
{
	obj_t* o;
	obj_t* o_1;
	obj_t* o_2;
}	tri_t;

//From data_utility.h
//The data structure for storing the statistics.
typedef struct coskq_stat
{
	//cost.
	B_KEY_TYPE	aver_cost;
    B_KEY_TYPE	aver_weight;
	B_KEY_TYPE	aver_B;
	
	//time.
	float		q_time;
	float		irtree_build_time;

	//memory.
	float		memory_v;
	float		memory_max;
	
	float		tree_memory_v;
	float		tree_memory_max;

	//float		data_memory_v;
	//float		data_memory_max;
	
	//algorithm related.
	float		n_1_sum;
	float		achi_sum;

	float		O_t_size_aver;
	float		O_t_size_sum;
	float		O_t_num;

	float		psi_n_aver;
	float		psi_n_sum;

	float		O_simp_size;

	//Cao-Exact specific.
	float		node_set_n;
	float		feasible_set_n;

	//Cao-Appro2
	float		n_k;

	//Additional.
//	float		ratio_min;
//	float		ratio_max;
	float		ratioa_aver;
	float		ratiob_aver;
	float		ratiob2_aver;
//	float		ratio_dev;

}	coskq_stat_t;

#endif
