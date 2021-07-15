//
//  cdcoskq.hpp
//  CoSKQ
//
//  Created by Harry on 23/4/2019.
//  Copyright © 2019 Harry. All rights reserved.
//

#ifndef cdcoskq_hpp
#define cdcoskq_hpp

#include <stdio.h>
#include "costenum.h"
//#include "cao_alg_new.h"
#include "bst.h"
#include "data_struct.h"
#include "b_heap.h"

#define EPSILON				0.0000001

extern B_KEY_TYPE MAX_DIA;
extern W_TYPE MAX_WEIGHT;

extern bst_t* IF_global;

extern int cost_tag;
extern int w_tag;

obj_set_t* Exact(query_t* q, int m_opt, double B);

obj_set_t* Appro(query_t* q, int m_opt, double B);

obj_set_t* ApproAdapt(query_t* q, int m_opt, double B);

obj_set_t* ExactBaseline(query_t* q, int m_opt, double B);

obj_set_t* RadIntMethods(query_t* q, const int alg_opt, const int m_opt, const double B);

obj_set_t* localOptimalFinding(int m_opt, B_KEY_TYPE B, loc_t* loc_q, psi_t* psi_v, obj_set_t* S_cur, obj_node_t* obj_node_v_t, W_TYPE curWeight, bool callFromPlus, W_TYPE start_x, W_TYPE end_x);

obj_set_t* localOptimalFinding_new(int m_opt, B_KEY_TYPE B, loc_t* loc_q, psi_t* psi_v, obj_set_t* S_cur, obj_node_t* obj_node_v_t, W_TYPE curWeight);

obj_set_t* baselineSearch(int m_opt, B_KEY_TYPE B, query_t* q, psi_t* psi_v, obj_set_t* S, obj_t* o_t, obj_set_t*& G_curBest, W_TYPE& curWeight);

bool search_sub(int m_opt, B_KEY_TYPE B, bst_t* IF_S, obj_set_t* G, query_t* q, B_KEY_TYPE d_o_q, W_TYPE weight, obj_set_t*& curBest, W_TYPE& curWeight, bool& foundResult);

bool search_sub2(int m_opt, B_KEY_TYPE B, bst_t* IF_S, obj_set_t* G, query_t* q, obj_set_t*& curBest, W_TYPE& curWeight);

obj_set_t* APPRO2(query_t* q, obj_set_t* S);

bst_t* build_Tab(psi_t* psi_v, obj_set_t* obj_set_v);

bst_t* build_Tab_min_weight(psi_t* psi_v, obj_set_t* obj_set_v, W_TYPE seed_weight);

bool get_min_weight_set(bst_node_t* bst_node_v, obj_set_t* obj_set_v);

bool check_Tab_not_zero(bst_node_t* bst_node_v);

bool check_Tab_all_one(bst_node_t* bst_node_v);

bool check_critical_obj_in_Tab(bst_node_t* bst_node_v, obj_t* obj_v);

b_heap_t* heap_sort_obj_set_weight( obj_set_t* obj_set_v);

bst_node_t* findMostInfreqKey(bst_t* IF_bst_v, psi_t* psi_v);

W_TYPE comp_weight(int weight_tag, obj_set_t* obj_set_v);

W_TYPE comp_weight_max(obj_set_t* obj_set_v);

W_TYPE comp_weight_sum(obj_set_t* obj_set_v);

void range_query_sub(node_t* node_v, disk_t* disk_v1, disk_t* disk_v2, obj_set_t*& obj_set_v, query_t* q);

obj_set_t* range_query(disk_t* disk_v1, disk_t* disk_v2, query_t* q);

void retrieve_sub_tree_within_weight( node_t* node_v, obj_set_t* &obj_set_v, query_t* q, W_TYPE min_weight);

void range_query_within_weight_sub(node_t* node_v, disk_t* disk_v1, disk_t* disk_v2, obj_set_t*& obj_set_v, query_t* q, W_TYPE max_weight);

obj_set_t* range_query_within_weight(disk_t* disk_v1, disk_t* disk_v2, query_t* q, W_TYPE max_weight);

void retrieve_sub_tree_ellipse_within_weight(node_t* node_v, obj_set_t*& obj_set_v, query_t* q,
											 W_TYPE max_weight);
void range_query_ellipse_within_weight_sub(node_t* node_v, disk_t* disk_v1, obj_set_t*& obj_set_v, query_t* q, W_TYPE max_weight, ellipse_t* ellipse_v);

obj_set_t* range_query_ellipse_within_weight(disk_t* disk_v1, query_t* q, W_TYPE max_weight, ellipse_t* ellipse_v);

BIT_TYPE has_key_node_within_weight( node_t* node_v, KEY_TYPE key, W_TYPE weight);

BIT_TYPE is_relevant_node_within_weight( node_t* node_v, query_t* q, W_TYPE weight);

bst_t* const_IF_sorted( obj_set_t* obj_set_v, psi_t* psi_v);

void insert_IF_obj(bst_t* IF_v, obj_t* obj_v);

double nthHarmonic(int N);

bool IF_obj_check(bst_node_t* bst_node_v);

W_TYPE findOptMaxWeightLB(bst_t* IF_bst_v, query_t* q, B_KEY_TYPE B);

void print_IF( bst_t* T);

void print_IF_node(bst_node_t* bst_node_v);

obj_set_t* comp_LB(query_t* q, disk_t* disk_q_B,B_KEY_TYPE& distLB);

obj_set_t* removeObjDominated_new(obj_set_t* O_t, psi_t* psi_v, bool returnRemoved);

obj_set_t** filterAndDistributeObj(obj_set_t* O_t, psi_t* psi_v, W_TYPE max_weight);

obj_set_t* localOptimalFinding_first(int m_opt, B_KEY_TYPE B, query_t* q, psi_t* psi_v, obj_set_t* S_cur, obj_t* obj_v,  bst_t*& IF_S, W_TYPE curWeight, bool& foundResult);

obj_set_t* localOptimalFinding_second(int m_opt, B_KEY_TYPE B, query_t* q, psi_t* psi_v, obj_set_t* S_cur, obj_t* obj_key, bst_t* IF_S, W_TYPE curWeight);

obj_set_t* localOptimalFinding_second_sum(int m_opt, B_KEY_TYPE B, query_t* q, psi_t* psi_v, obj_set_t* S_cur, obj_t* obj_key, bst_t* IF_S, W_TYPE curWeight, bool& foundResult);

bool is_enclosed( range* MBR, ellipse_t* ellipse_v);

B_KEY_TYPE calc_maxDist(range* MBR, loc_t* loc_v1, loc_t* loc_v2);

W_TYPE findOptSumWeightLB(bst_t* IF_bst_v, query_t* q, B_KEY_TYPE B);

void update_Tab(bst_node_t* bst_node_v, obj_t* obj_v);

obj_set_t* localOptimalFinding_Max(double B, obj_set_t *&S_cur, int curWeight, int m_opt, obj_t *obj_v, query_t *q, query_t *q_new);

obj_set_t* locaOptimalFinding_Sum(double B, obj_set_t *&S_cur, int curWeight, int m_opt, obj_t *obj_v, query_t *q, query_t *q_new);

obj_set_t* pickDominateObj(obj_set_t* O_t, psi_t* psi_v);

bool APPRO1_addver2(query_t* q, obj_set_t* S, obj_set_t* G, bst_t* Tab);

obj_set_t* localFeasibleFinding_Adapt(obj_t* o, query_t* q, disk_t* disk_v);

obj_set_t* localFeasibleFinding_Max(obj_set_t *&S_cur, obj_t *obj_v, query_t *q_new);

obj_set_t* findMinWeightSet(bst_t* IF_bst_v, query_t* q, B_KEY_TYPE B);

int is_relevant_obj(obj_t* obj_v, psi_t* psi_v);

void release_IF(bst_t* T);

void release_IF_sub(bst_node_t* x);

bool check_dist_constraint(obj_set_t* obj_set_v, obj_t* obj_v, B_KEY_TYPE d);

bool check_pairwise_dist_constraint(obj_set_t* obj_set_v, B_KEY_TYPE d);

bool isDominated(obj_t* obj_v, obj_t* obj_v2, psi_t* psi_v);

int number_intersection(k_node_t* k_head1, k_node_t* k_head2);

psi_t* psi_exclusion(psi_t* psi_v1, obj_t* obj_v);

void psi_insert(psi_t* psi_v, k_node_t* k_head);

#endif /* cdcoskq_hpp */
