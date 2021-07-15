//
//  cdcoskq.cpp
//  CoSKQ
//
//  Created by Harry on 23/4/2019.
//  Copyright © 2019 Harry. All rights reserved.
//

#include "cdcoskq.h"

obj_set_t* Exact(query_t* q, int m_opt, double B)
{
    return RadIntMethods(q, 1, m_opt, B);
}

obj_set_t* Appro(query_t* q, int m_opt, double B)
{
    return RadIntMethods(q, 2, m_opt, B);
}

obj_set_t* ExactBaseline(query_t* q, int m_opt, double B)
{
    return RadIntMethods(q, 4, m_opt, B);
}

/*
 m_opt
 =1 : Cao-A2
 =2 : Long-A
 */
obj_set_t* ApproAdapt(query_t* q, int m_opt, double B)
{
    //    bst_t* IF_v;
    obj_set_t *G_opt, *O, *G;
    //    obj_node_t* obj_node_v;
    disk_t *disk_v1, *disk_v2, *disk_q_B;
    query_t* q_new;
    loc_t* loc_v;
    W_TYPE curWeight;
    double d_o_q;
    k_node_t* k_head_1;
    obj_t* obj_v;
    B_KEY_TYPE LB, curCost;
    //    bool foundResult;
    int i = 0;
    W_TYPE optWeightLB;

    b_heap_t* R_heap;
    int top;
    obj_t *o, *o_next;

    G_opt = NULL;
    curWeight = w_tag == 1 ? MAX_WEIGHT : MAX_WEIGHT * q->psi_v->key_n;
    curCost = INFINITY;

    // Construct a disk center at q with radius B
    disk_q_B = alloc_disk(2);
    set_disk(disk_q_B, q->loc_v, B);
    G_opt = comp_LB(q, disk_q_B, LB); //note that this G_opt may not satisfy distance constraint

    //update curWeight by N(q).
    if (G_opt != NULL) {
        curCost = comp_cost(cost_tag, G_opt, q);
        if (curCost <= 1 * B) {
            curWeight = comp_weight(w_tag, G_opt);
        }
    }

    //Perform range query and construct inverted lists
    //	obj_set_B = range_query(disk_q_B, q);
    //	obj_set_B = range_query_within_weight(disk_q_B, NULL,q, curWeight);
    //	IF_v = const_IF(obj_set_B, q->psi_v);

    if (w_tag == 1) {
        obj_set_t* G2 = findMinWeightSet(IF_global, q, B);
        if (G2 != NULL) {
            W_TYPE G2weight = comp_weight(w_tag, G2);
            if (G2weight <= curWeight && comp_cost(cost_tag, G2, q) <= B) {
                //				printf("-here\n");
                //				release_obj_set(obj_set_B);
                release_obj_set(G_opt);
                //				release_IF(IF_v);
                G_opt = G2;
                curWeight = G2weight;
                return G_opt;
            }
            release_obj_set(G2);
        }
    }

    //---
    if (w_tag == 1) {
        optWeightLB = findOptMaxWeightLB(IF_global, q, B);
    } else //w_tag==2
        optWeightLB = findOptSumWeightLB(IF_global, q, B);

    //---
    //    printf("optWeightLB:%d\tcurWeight:%d\n", optWeightLB, curWeight);

    // m_opt==1: (Radius Mehod) for each obj with t_inf
    // m_op2==2: (Intersection Method) for each object in D(q,B)
    if (m_opt == 1) {
        O = copy_obj_set(findMostInfreqKey(IF_global, q->psi_v)->p_list_obj);
    } else // m_opt ==2
    {
        //		O = range_query_within_weight(disk_q_B, NULL, q, curWeight);
        O = range_query(disk_q_B, q);

        disk_t* disk_LB = alloc_disk(2);
        set_disk(disk_LB, q->loc_v, LB);

        ///exclude the objects in region_u that are inside disk_l
        obj_exclusion_disk(O, disk_LB);
        release_disk(disk_LB);
    }
    //    printf("O->obj_n:%d\n", O->obj_n);

    //Sort the objects in R by their distance in ascending order
    R_heap = heap_sort_obj_set(O, q);
    //	R_heap = heap_sort_obj_set_weight(O);

    top = b_h_get_top(R_heap);

    o_next = R_heap->obj_arr[top].obj_v;

    while (true) {
        if (o_next == NULL)
            break;

        o = o_next;

        top = b_h_get_top(R_heap);
        if (top == 0)
            o_next = NULL;
        else
            o_next = R_heap->obj_arr[top].obj_v;

        obj_v = o;
        //	printf("i:%d\t obj_node_v:%d\n", i++, obj_node_v->obj_v->id);
        // Pruning: skip this object if its weight >= w(G_curBest).
        if (obj_v->weight >= curWeight) {
			continue;
        }

        loc_v = get_obj_loc(obj_v);
        d_o_q = calc_dist_loc(loc_v, q->loc_v);

        // Pruning: skip this object if d_o_q > B or d_o_q <LB
        //objs in m_opt ==2 are already in D(q,B)
        if ((m_opt == 1 && d_o_q > B) || (m_opt == 2 && d_o_q < LB)) {
            release_loc(loc_v);
            continue;
        }

        if (m_opt == 2 && d_o_q > curCost) {
            release_loc(loc_v);
            break;
        }


        // Construct a new query instance q_new.
        // with loc = o_t.loc and \psi = q.\psi - o.\psi
        q_new = alloc_query();
        q_new->loc_v = loc_v;
        k_head_1 = key_exclusion(q->psi_v->k_head, obj_v->k_head);
        q_new->psi_v = const_psi(k_head_1);

        disk_v2 = NULL;
        if (m_opt == 2) {
            //disk_v2 == Disk(q,d(o,q0))
            disk_v2 = alloc_disk(2);
            set_disk(disk_v2, q->loc_v, d_o_q);
        }

        // Part 2.

        //m_opt is indicated by disk_v2
        G = localFeasibleFinding_Adapt(obj_v, q_new, disk_v2);

        if (G != NULL) {
            add_obj_set_entry(obj_v, G);
            //update current best solution
            W_TYPE G_weight = comp_weight(w_tag, G);
            //					printf("G_weight:%d\tcurWeight:%d\toptWeightLB:%d\t LB:%lf\n",G_weight,curWeight,optWeightLB, optWeightLB*nthHarmonic(q->psi_v->key_n));
            B_KEY_TYPE G_cost = comp_cost(cost_tag, G, q);
            if (G_cost < curCost) {
                curCost = G_cost;
            }
            if (G_cost <= 1 * B && G_weight <= curWeight) {
                release_obj_set(G_opt);
                G_opt = G;
                curWeight = G_weight;
                if (curWeight == optWeightLB)
                    break;
            }
        }

        release_disk(disk_v2);
        release_query(q_new);

    } //end while

    release_disk(disk_q_B);
    release_b_heap(R_heap);
    release_obj_set(O);

    //	release_IF(IF_v);

    return G_opt;
}

/*
 * Find NN of o for each key in q
 *
 * disk_v == NULL indicate the whole space -> Cao-Appro
 * otherwise disk_v = Disk(q,d(o,q)) -> Long-Appro
 */
obj_set_t* localFeasibleFinding_Adapt(obj_t* o, query_t* q, disk_t* disk_v)
{
    obj_set_t* S;
    obj_t* obj_v;
    k_node_t* k_node_v;
    loc_t* loc_v;

    S = alloc_obj_set();
    loc_v = get_obj_loc(o);

    k_node_v = q->psi_v->k_head->next;
    while (k_node_v != NULL) {
        obj_v = const_NN_key(loc_v, k_node_v->key, disk_v);

        if (obj_v == NULL) {
            printf("No solution exists!\n");
            release_obj_set(S);
            release_loc(loc_v);
            return NULL;
        }

        add_obj_set_entry(obj_v, S);

        k_node_v = k_node_v->next;
    }
    release_loc(loc_v);

    return S;
}

obj_set_t* localFeasibleFinding_Max(obj_set_t*& S_cur, obj_t* obj_v, query_t* q_new)
{
    obj_set_t* G;
    obj_set_t** obj_set_set;
    bst_t* Tab;

    //					print_obj_set(S_cur, stdout);
    //					printf("S_cur->obj_n:%d\n", S_cur->obj_n);

    //S_cur is updated in this function
    obj_set_set = filterAndDistributeObj(S_cur, q_new->psi_v, obj_v->weight);

    //----
    Tab = build_Tab(q_new->psi_v, S_cur);
    if (check_Tab_not_zero(Tab->root)) {
        bst_release(Tab);
        G = S_cur;
        S_cur = NULL;
        for (int x = 0; x < 10; x++)
            release_obj_set(obj_set_set[x]);
        free(obj_set_set);
        return G;
    }
    //----
    G = S_cur;
    S_cur = NULL;
    bool feasible = false;
    int j = 0;
    while (!feasible && j < 10) {

        release_obj_set(S_cur); //IF_S already stored all objects in S_cur in each iteration
        S_cur = obj_set_set[j];
        obj_set_set[j] = NULL;
        //						printf("j:%d\tS_cur:%d\n",j,S_cur->obj_n);
        //						print_obj_set(S_cur, stdout);
        if (S_cur->obj_n > 0)
            feasible = APPRO1_addver2(q_new, S_cur, G, Tab);
        j++;
    }

    bst_release(Tab);
    for (; j < 10; j++)
        release_obj_set(obj_set_set[j]);
    free(obj_set_set);

    return G;
}

obj_set_t* localOptimalFinding_Max(double B, obj_set_t*& S_cur, int curWeight, int m_opt, obj_t* obj_v, query_t* q, query_t* q_new)
{
    obj_set_t* G;
    obj_set_t** obj_set_set;
    bst_t* IF_S;
    bool foundResult;

    //					print_obj_set(S_cur, stdout);
    //					printf("S_cur->obj_n:%d\n", S_cur->obj_n);

    foundResult = false;

    //S_cur is updated in this function
    obj_set_set = filterAndDistributeObj(S_cur, q_new->psi_v, obj_v->weight);
    G = localOptimalFinding_first(m_opt, B, q, q_new->psi_v, S_cur, obj_v, IF_S, curWeight, foundResult);

    //process a bucket in each iteration until (optimal) solution is found
    int j = 0;
    while (G == NULL && j < 10) {

        release_obj_set(S_cur); //IF_S already stored all objects in S_cur in each iteration
        S_cur = obj_set_set[j];
        obj_set_set[j] = NULL;
        //						printf("j:%d\tS_cur:%d\n",j,S_cur->obj_n);
        //						print_obj_set(S_cur, stdout);
        if (S_cur != NULL && S_cur->obj_n > 0)
            G = localOptimalFinding_second(m_opt, B, q, q_new->psi_v, S_cur, obj_v, IF_S, curWeight);
        j++;
    }
    release_IF(IF_S);
    for (; j < 10; j++)
        release_obj_set(obj_set_set[j]);
    free(obj_set_set);

    return G;
}

obj_set_t* locaOptimalFinding_Sum(double B, obj_set_t*& S_cur, int curWeight, int m_opt, obj_t* obj_v, query_t* q, query_t* q_new)
{
    obj_set_t* G;
    obj_set_t* S_cur2;
    bst_t* IF_S;
    bool foundResult;
    //					printf("i:%d\t\toriginal S_cur:%d\n",i,S_cur->obj_n);

    //flag to indicate local optimal is found
    //if true mean do not need to further process this iteration
    //since G could be NULL if weight of local opt found is larger than curWeight
    //so check G == NULL is not sufficient
    foundResult = false;

    S_cur2 = removeObjDominated_new(S_cur, q_new->psi_v, true);
    G = localOptimalFinding_first(m_opt, B, q, q_new->psi_v, S_cur, obj_v, IF_S, curWeight, foundResult);

    release_obj_set(S_cur);
    S_cur = S_cur2;

    while (G == NULL && S_cur->obj_n > 0 && !foundResult) {
        S_cur2 = removeObjDominated_new(S_cur, q_new->psi_v, true);

        //printf("S_cur:%d\tS_cur2:%d\tc:%d\n", S_cur->obj_n, S_cur2->obj_n, c);

        if (S_cur->obj_n > 0)
            G = localOptimalFinding_second_sum(m_opt, B, q, q_new->psi_v, S_cur, obj_v, IF_S, curWeight, foundResult);
        release_obj_set(S_cur);
        S_cur = S_cur2;
    }

    release_IF(IF_S);

    return G;
}

/*
 * alg_opt:
 *	1: exact
 *	2: appro
 *	4: baseline exact
 *
 * m_opt:
 *	1: radius method
 *	2: intersection method
 *
 */
obj_set_t* RadIntMethods(query_t* q, const int alg_opt, const int m_opt, const double B)
{

    //    bst_t* IF_v;

    obj_set_t *G_opt, *S_cur, *O, *G, *G2, *obj_set_infreq;
    obj_node_t* obj_node_infreq;
    //    obj_node_t* obj_node_v;
    disk_t *disk_v1, *disk_v2, *disk_q_B;
    query_t* q_new;
    loc_t* loc_v;
    W_TYPE curWeight;
    double d_o_q;
    k_node_t* k_head_1;
    obj_t* obj_v;
    B_KEY_TYPE LB;
    //    bool foundResult;
    int i = 0;
    W_TYPE optWeightLB;

    b_heap_t* R_heap;
    int top;
    obj_t *o, *o_next;

    G_opt = NULL;
    curWeight = w_tag == 1 ? MAX_WEIGHT : MAX_WEIGHT * q->psi_v->key_n;

    // Construct a disk center at q with radius B
    disk_q_B = alloc_disk(2);
    set_disk(disk_q_B, q->loc_v, B);
    G_opt = comp_LB(q, disk_q_B, LB); //NN set, note that this G_opt may not satisfy distance constraint

    //update curWeight by N(q).
    if (G_opt != NULL) {
        if (alg_opt == 1) { //exact
            if (comp_cost(cost_tag, G_opt, q) <= B)
                curWeight = comp_weight(w_tag, G_opt);
        } else if (alg_opt == 2) { //approx
            if (m_opt == 1) {
                if ((cost_tag == 1 && comp_cost(cost_tag, G_opt, q) <= 1.5 * B)) {
                    curWeight = comp_weight(w_tag, G_opt);
                } else if ((cost_tag == 2 && comp_cost(cost_tag, G_opt, q) <= 1.5 * B)) {
                    curWeight = comp_weight(w_tag, G_opt);
                }
            } else //m_opt==2
            {
                if (cost_tag == 1 && comp_cost(cost_tag, G_opt, q) <= 1.375 * B)
                    curWeight = comp_weight(w_tag, G_opt);
                else if (cost_tag == 2 && comp_cost(cost_tag, G_opt, q) <= sqrt(3) * B)
                    curWeight = comp_weight(w_tag, G_opt);
            }
        }
    }

    //Pruning based on Min. Weight set
    G2 = findMinWeightSet(IF_global, q, B);
    if (G2 != NULL) {
        W_TYPE G2weight = comp_weight(w_tag, G2);
        if (G2weight <= curWeight) {
            if (alg_opt == 1) { //exact
                if (comp_cost(cost_tag, G2, q) <= B) {
                    if (w_tag == 1) {
                        //                        printf("---return here\n");
                        release_obj_set(G_opt);
                        return G2;
                    } else {
                        curWeight = G2weight;
                        release_obj_set(G_opt);
                        G_opt = G2;
                        G2 = NULL;
                    }
                }
            } else if (alg_opt == 2) { //approx
                if (m_opt == 1) {
                    if ((cost_tag == 1 && comp_cost(cost_tag, G2, q) <= 1.5 * B)
                        || (cost_tag == 2 && comp_cost(cost_tag, G2, q) <= 1.5 * B)) {
                        if (w_tag == 1) {
                            //                                                        printf("--return here\n");
                            release_disk(disk_q_B);
                            release_obj_set(G_opt);
                            //                            release_IF(IF_v);

                            return G2;
                        } else {
                            curWeight = G2weight;
                            release_obj_set(G_opt);
                            G_opt = G2;
                            G2 = NULL;
                        }
                    }
                } else //m_opt==2
                {
                    if ((cost_tag == 1 && comp_cost(cost_tag, G2, q) <= 1.375 * B)
                        || (cost_tag == 2 && comp_cost(cost_tag, G2, q) <= sqrt(3) * B)) {
                        if (w_tag == 1) {
                            //                            printf("-return here\n");
                            release_disk(disk_q_B);
                            release_obj_set(G_opt);
                            //                            release_IF(IF_v);
                            return G2;
                        } else {
                            curWeight = G2weight;
                            release_obj_set(G_opt);
                            G_opt = G2;
                            G2 = NULL;
                        }
                    }
                }
            }
        }
        release_obj_set(G2);
    } else {
        printf("G2 not found \n");
        exit(-1);
    }

    //---

        //Pruning based on Greedy Set Cover on whole circle
        if ( w_tag == 2) {
			obj_set_t* obj_set_B = range_query_within_weight(disk_q_B, NULL, q, curWeight);
			
            G2 = APPRO2(q, obj_set_B);
            W_TYPE G2weight = comp_weight(w_tag, G2);
            //        		printf("optWeightLB2:%d\n",G2weight);
            if (alg_opt == 1) {
                if (comp_cost(cost_tag, G2, q) <= B) {
    
                    //exact, update curWeight
                    if (G2weight <= curWeight) {
                        curWeight = G2weight;
                        release_obj_set(G_opt);
                        G_opt = G2;
                        G2 = NULL;
                    }
                }
            } else if (alg_opt == 2) {
                //approx, return here
                if (m_opt == 1) {
                    if ((cost_tag == 1 && comp_cost(cost_tag, G2, q) <= 1.5 * B)
                        || (cost_tag == 2 && comp_cost(cost_tag, G2, q) <= 1.5 * B)) {
    
                        //				printf("-return\n");
                        release_obj_set(G_opt);
						release_obj_set(obj_set_B);
                        return G2;
                    }
                } else {//m_opt==2
                    if ((cost_tag == 1 && comp_cost(cost_tag, G2, q) <= 1.375 * B)
                        || (cost_tag == 2 && comp_cost(cost_tag, G2, q) <= sqrt(3) * B)) {
    
                        //				printf("-return\n");
                        release_obj_set(G_opt);
						release_obj_set(obj_set_B);
                        return G2;
                    }
                }
            }
			
            release_obj_set(G2);
			release_obj_set(obj_set_B);
        }

    //---
    if (w_tag == 1)
        optWeightLB = findOptMaxWeightLB(IF_global, q, B);
    else { //w_tag==2
        optWeightLB = findOptSumWeightLB(IF_global, q, B);
        //        printf("optWeightLB1:%d\n", optWeightLB);
    }
    //---

    //    printf("optWeightLB:%d\tcurWeight:%d\n", optWeightLB, curWeight);

    if (curWeight == optWeightLB) {
        release_disk(disk_q_B);

        //        release_IF(IF_v);
        return G_opt;
    }

    // m_opt==1: (Radius Mehod) for each obj with t_inf
    // m_op2==2: (Intersection Method) for each object in D(q,B)
    if (m_opt == 1) {
        O = copy_obj_set(findMostInfreqKey(IF_global, q->psi_v)->p_list_obj);

        //		O = copy_obj_set(find_largest_min_dist_key(IF_global,q)->p_list_obj);

    } else // m_opt ==2
    {
        O = range_query_within_weight(disk_q_B, NULL, q, curWeight);

        disk_t* disk_LB = alloc_disk(2);
        set_disk(disk_LB, q->loc_v, LB);

        ///exclude the objects in region_u that are inside disk_l
        obj_exclusion_disk(O, disk_LB);
        release_disk(disk_LB);
    }
    //    printf("O->obj_n:%d\n", O->obj_n);

    //Sort the objects in R by their weight in ascending order
    R_heap = heap_sort_obj_set_weight(O);
    top = b_h_get_top(R_heap);

    o_next = R_heap->obj_arr[top].obj_v;

    obj_set_infreq = findMostInfreqKey(IF_global, q->psi_v)->p_list_obj; //no memory release is needed

    while (true) {
        if (o_next == NULL)
            break;

        o = o_next;

        top = b_h_get_top(R_heap);
        if (top == 0)
            o_next = NULL;
        else
            o_next = R_heap->obj_arr[top].obj_v;

        obj_v = o;
        //	printf("i:%d\t obj_node_v:%d\n", i++, obj_node_v->obj_v->id);
        // Pruning: skip this object if its weight >= w(G_curBest).
        if (obj_v->weight >= curWeight) {
            break;
            //            continue;
        }

        loc_v = get_obj_loc(obj_v);
        d_o_q = calc_dist_loc(loc_v, q->loc_v);

        // Pruning: skip this object if d_o_q < LB
        if ((m_opt == 1 && d_o_q > B) || (m_opt == 2 && d_o_q < LB)) {
            release_loc(loc_v);
            continue;
        }

        //---
        //skip obj if its S_cur do not have obj cover most infrequent keyword
        //work (i.e., reduced number of objects to be processed)
        //
        if (m_opt == 2) {
            bool feasible = false;
            obj_node_infreq = obj_set_infreq->head->next;
            while (obj_node_infreq != NULL) {
                loc_t* loc_infreq = get_obj_loc(obj_node_infreq->obj_v);
                //check whether obj_node_infreq located in the intersection of two disks
                if (calc_dist_loc(loc_v, loc_infreq) <= B - d_o_q
                    && calc_dist_loc(q->loc_v, loc_infreq) <= d_o_q) {
                    feasible = true;
                    break;
                }
                release_loc(loc_infreq);
                obj_node_infreq = obj_node_infreq->next;
            }
            if (!feasible) {

                continue;
            }
        }

        //        printft("doq:%lf\tB:%lf\t%lf\tw:%d\n", d_o_q, B, B - d_o_q, obj_v->weight);

        //---

        

        // Construct a new query instance q_new.
        // with loc = o_t.loc and \psi = q.\psi - o.\psi
        q_new = alloc_query();
        q_new->loc_v = loc_v;
        k_head_1 = key_exclusion(q->psi_v->k_head, obj_v->k_head);
        q_new->psi_v = const_psi(k_head_1);

        // Find the set S_cur
        disk_v1 = alloc_disk(2);
        if (cost_tag == 1) //MaxSum
            set_disk(disk_v1, loc_v, B - d_o_q);
        else //cost_tag==2 //Dia
            set_disk(disk_v1, loc_v, B);

        if (m_opt == 1) //radius
        {
            if (alg_opt == 4) //no ellipse pruning and weight pruning in baseline exact
                S_cur = range_query(disk_v1, disk_q_B, q_new);
            else { //alg_opt ==1 || alg_opt==2
                if (cost_tag == 1) {
                    ellipse_t* ellipse_v = const_ellipse(loc_v, q->loc_v, B);
                    S_cur = range_query_ellipse_within_weight(disk_v1, q_new, curWeight, ellipse_v);
                    release_ellipse(ellipse_v);
                } else { //cost_tag==2
                    S_cur = range_query_within_weight(disk_v1, disk_q_B, q_new, curWeight);
                }
            }
        } else // m_opt==2 intersection
        {
            disk_v2 = alloc_disk(2);
            set_disk(disk_v2, q->loc_v, d_o_q);
            //			S_cur = find_intersest_obj(disk_v1, disk_v2, obj_set_B, q_new, curWeight);
            S_cur = range_query_within_weight(disk_v1, disk_v2, q_new, curWeight);
			release_disk(disk_v2);
        }
		stat_v.O_t_size_sum+= S_cur->obj_n;
		stat_v.O_t_num++;
		
        // Part 2.
        if (alg_opt == 1) // exact
        {
            if (FeasibilityCheck(S_cur, q_new->psi_v)) {
                if (w_tag == 1) //Max Weight
                    G = localOptimalFinding_Max(B, S_cur, curWeight, m_opt, obj_v, q, q_new);
                else //w_tag==2
                    G = locaOptimalFinding_Sum(B, S_cur, curWeight, m_opt, obj_v, q, q_new);

                //update current best solution
                if (G != NULL) {
                    W_TYPE G_weight = comp_weight(w_tag, G);
                    //					printf("G_weight:%d\tcurWeight:%d\toptWeightLB:%d\n",G_weight,curWeight,optWeightLB);

                    if (G_weight <= curWeight) {
                        release_obj_set(G_opt);
                        G_opt = G;
                        curWeight = G_weight;
                        if (curWeight == optWeightLB) {
                            //optimal found, terminate here
                            if (S_cur != NULL)
                                release_obj_set(S_cur);
                            release_query(q_new);
                            release_disk(disk_v1);
                            break;
                        }
                    }
                }
            }
        } else if (alg_opt == 2) // approx
        {
            //m_opt==2 && d_o_q<=B/3, S_cur must include N(q) and must be feasible
            if ((m_opt == 2 && d_o_q <= B / 3) || (w_tag==1|| FeasibilityCheck(S_cur, q_new->psi_v))) {
                if (w_tag == 1) {
					//check feasibility is done inside
                    G = localFeasibleFinding_Max(S_cur, obj_v, q_new);

                } else // weight_tag ==2
                {

                    //speed up by only consider dominate objects
                    //					printf("S_cur->obj_n:%d\n",S_cur->obj_n);
                    obj_set_t* dominateObjs = pickDominateObj(S_cur, q_new->psi_v);
                    G = APPRO2(q_new, dominateObjs);
                    release_obj_set(dominateObjs);
                }

                if (G != NULL) {
                    add_obj_set_entry(obj_v, G);
                    //update current best solution
                    W_TYPE G_weight = comp_weight(w_tag, G);
                    //					printf("G_weight:%d\tcurWeight:%d\toptWeightLB:%d\t LB:%lf\n",G_weight,curWeight,optWeightLB, optWeightLB*nthHarmonic(q->psi_v->key_n));

                    if (G_weight <= curWeight) {
                        release_obj_set(G_opt);
                        G_opt = G;
                        curWeight = G_weight;
                        //w_tag==2 early stopping checking:
                        //-1 here to have better approximation in practice
                        //theorically no need
                        if (curWeight == optWeightLB || (w_tag == 2 && curWeight <= optWeightLB * (nthHarmonic(q->psi_v->key_n) - 1))) {
                            //optimal found, terminate here
                            if (S_cur != NULL)
                                release_obj_set(S_cur);
                            release_query(q_new);
                            release_disk(disk_v1);
                            break;
                        }
                    }
                }
            }
        } else if (alg_opt == 4) // baselineexact
        {
            if (FeasibilityCheck(S_cur, q_new->psi_v))
                baselineSearch(m_opt, B, q, q_new->psi_v, S_cur, obj_v, G_opt, curWeight);
        }

    E:
        if (S_cur != NULL)
            release_obj_set(S_cur);
        release_query(q_new);
        release_disk(disk_v1);

        //        obj_node_v = obj_node_v->next;
    } //end while

    release_disk(disk_q_B);
    release_b_heap(R_heap);
    release_obj_set(O);

    //    release_IF(IF_v);

    return G_opt;
}

/*
 *	Find an exact solution (i.e., a set with min weight, cost<B, cover all keywords) if exists
 *	return NULL otherwise
 *
 *	@o_t is the object must be included
 *	@S_cur is the set of objects to be enumerated
 *	@q is the original query (for retreive q->loc_v)
 *	@psi_v = q.\psi - o_t.\psi
 *
 * first: build IF + initial search
 **/
obj_set_t* localOptimalFinding_first(int m_opt, B_KEY_TYPE B, query_t* q, psi_t* psi_v, obj_set_t* S_cur, obj_t* obj_v, bst_t*& IF_S, W_TYPE curWeight, bool& foundResult)
{
    obj_set_t *G, *G_o;
    loc_t* loc_v;
    B_KEY_TYPE d_o_q;

    G_o = NULL;

    //Construct "initial" IF_S with S_cur
    IF_S = const_IF(S_cur, psi_v);

    if (IF_obj_check(IF_S->root)) {
        G = alloc_obj_set();
        add_obj_set_entry(obj_v, G);

        loc_v = get_obj_loc(obj_v);
        d_o_q = calc_dist_loc(loc_v, q->loc_v);
        release_loc(loc_v);

        //check whether this IF_S can find a solution
        search_sub(m_opt, B, IF_S, G, q, d_o_q, obj_v->weight, G_o, curWeight, foundResult);
        release_obj_set(G);
    }
    return G_o;
}

/*
 * iterate each object o in S_cur
 * 1. check if o + objects in IF_S can cover psi and feasible and within cost
 * 2. if not, add to IF_S
 * 4. ALL objects in S_cur need to be iterated to find opt sol
 * 5. return null if no sol is found
 */
obj_set_t* localOptimalFinding_second(int m_opt, B_KEY_TYPE B, query_t* q, psi_t* psi_v, obj_set_t* S_cur, obj_t* obj_key, bst_t* IF_S, W_TYPE curWeight)
{
    obj_set_t *G, *G_o;
    B_KEY_TYPE d_o_q, d_o_q_new;
    bst_node_list_t* bst_node_list_v;
    W_TYPE weight;
    bool tag, foundResult;
    loc_t* loc_v;

    b_heap_t* R_heap;
    int top;
    obj_t *o, *o_next;

    G_o = NULL;
    G = alloc_obj_set();
    add_obj_set_entry(obj_key, G);

    loc_v = get_obj_loc(obj_key);
    d_o_q = calc_dist_loc(loc_v, q->loc_v);
    release_loc(loc_v);

    //Sort the objects in R by their weight.
    R_heap = heap_sort_obj_set_weight(S_cur);
    top = b_h_get_top(R_heap);

    o_next = R_heap->obj_arr[top].obj_v;
    while (true) {
        if (o_next == NULL)
            break;

        o = o_next;

        top = b_h_get_top(R_heap);
        if (top == 0)
            o_next = NULL;
        else
            o_next = R_heap->obj_arr[top].obj_v;

        // Update the IF_S
        bst_node_list_v = update_IF_obj(IF_S, o);

        //if IF_S contain "enough" object (i.e., all objects in IF_S is a feasible set)
        if (IF_obj_check(IF_S->root)) {

            // Update the S_0.
            // obj_v is added at the first place of G
            add_obj_set_entry(o, G);

            weight = comp_weight(w_tag, G);

            //update the cost.
            if (m_opt == 1) //Radius
            {
                d_o_q_new = calc_minDist(o->MBR, q->loc_v);
                d_o_q_new = d_o_q_new > d_o_q ? d_o_q_new : d_o_q;
            } else //intersection do not need to check
                d_o_q_new = d_o_q;

            tag = search_sub(m_opt, B, IF_S, G, q, d_o_q_new, weight, G_o, curWeight, foundResult);

            // Restore the S_0.
            remove_obj_set_entry(G);

            // Restore the IF_S.
            restore_IF_bst_node_list(IF_S, bst_node_list_v);
            release_bst_node_list(bst_node_list_v);

            if (w_tag == 1 && tag)
                break;

        } else {
            restore_IF_bst_node_list(IF_S, bst_node_list_v);
            release_bst_node_list(bst_node_list_v);
        }
        insert_IF_obj(IF_S, o);
    }
    //        printf("\n");
    //----
    release_obj_set(G);
    release_b_heap(R_heap);
    return G_o;
}

/*
 * iterate each object o in S_cur
 * 1. check if o + objects in IF_S can cover psi and feasible and within cost
 * 2. if not, add to IF_S
 * 4. ALL objects in S_cur need to be iterated to find opt sol
 * 5. return null if no sol is found
 */
obj_set_t* localOptimalFinding_second_sum(int m_opt, B_KEY_TYPE B, query_t* q, psi_t* psi_v, obj_set_t* S_cur, obj_t* obj_key, bst_t* IF_S, W_TYPE curWeight, bool& foundResult)
{
    obj_set_t *G, *G_o;
    B_KEY_TYPE d_o_q, d_o_q_new;
    obj_node_t* obj_node_v;
    bst_node_list_t* bst_node_list_v;
    W_TYPE weight;
    bool tag;
    loc_t* loc_v;

    b_heap_t* R_heap;
    int top;
    obj_t *o, *o_next;

    G_o = NULL;
    G = alloc_obj_set();
    add_obj_set_entry(obj_key, G);

    loc_v = get_obj_loc(obj_key);
    d_o_q = calc_dist_loc(loc_v, q->loc_v);
    release_loc(loc_v);


    obj_node_v = S_cur->head->next;
    while (obj_node_v != NULL) {
        o = obj_node_v->obj_v;

        // Update the IF_S.
        bst_node_list_v = update_IF_obj(IF_S, o);

        if (IF_obj_check(IF_S->root)) {

            // Update the S_0.
            // obj_v is added at the first place of S_0.
            add_obj_set_entry(o, G);

            weight = comp_weight(w_tag, G);

            //update the cost.
            if (m_opt == 1) //Radius
            {
                d_o_q_new = calc_minDist(o->MBR, q->loc_v);
                d_o_q_new = d_o_q_new > d_o_q ? d_o_q_new : d_o_q;
            } else //intersection not need to check
                d_o_q_new = d_o_q;

            tag = search_sub(m_opt, B, IF_S, G, q, d_o_q_new, weight, G_o, curWeight, foundResult);

            // Restore the S_0.
            remove_obj_set_entry(G);

            // Restore the IF_S.
            restore_IF_bst_node_list(IF_S, bst_node_list_v);
            release_bst_node_list(bst_node_list_v);

            if (w_tag == 1 && tag)
                break;

        } else {
            restore_IF_bst_node_list(IF_S, bst_node_list_v);
            release_bst_node_list(bst_node_list_v);
        }
        insert_IF_obj(IF_S, o);
        obj_node_v = obj_node_v->next;
    }
    //        printf("\n");
    //----
    release_obj_set(G);
    //	release_b_heap(R_heap);
    return G_o;
}

/*
 *	The sub-procedure of function "EnumerateSubset".
 *  curBest and curCost store the current best set and cost, respectively.
 *
 *  weight = the weight of the set to be found (for w1)
 *		= the LB of the weight of the set to be found (for w2)
 */
bool search_sub(int m_opt, B_KEY_TYPE B, bst_t* IF_S, obj_set_t* G, query_t* q, B_KEY_TYPE d_o_q, W_TYPE weight, obj_set_t*& curBest, W_TYPE& curWeight, bool& foundResult)
{
    obj_t* obj_v;
    bst_node_t* bst_node_v;
    obj_node_t* obj_node_v;
    bst_node_list_t* bst_node_list_v;
    W_TYPE weight_new;
    B_KEY_TYPE d, d_o_q_new;
    bool tag;

    if (IF_S->node_n == 0) {
        foundResult = true;
        //		printf("node_n==0\n");
        //update current best
        //        weight = comp_weight(w_tag, G);
        if (weight <= curWeight) {
            release_obj_set(curBest);
            curBest = copy_obj_set(G);
            curWeight = weight;
        }
        return true;
    }

    bst_node_v = IF_S->root;
    obj_node_v = bst_node_v->p_list_obj->head->next;
    while (obj_node_v != NULL) { // Pick an object.
        obj_v = obj_node_v->obj_v;
        //-----------------------------
        // Budget checking.
        if (m_opt == 1) {
            d_o_q_new = calc_minDist(obj_v->MBR, q->loc_v);
            if (d_o_q_new > d_o_q) {
                if (cost_tag == 1) {
                    //perform additional checking if d_o_q is updated
                    //since pairwise dist may violate the new (smaller) d
                    d = B - d_o_q_new;
                    if (!check_pairwise_dist_constraint(G, d)) {
                        obj_node_v = obj_node_v->next;
                        continue;
                    }
                }
            } else
                //keep the original d_o_q
                d_o_q_new = d_o_q;

        } else //m_opt==2
        {
            d_o_q_new = d_o_q;
        }
        // check if o can be added to G
        d = cost_tag == 1 ? B - d_o_q_new : B;
        if (!check_dist_constraint(G, obj_v, d)) {
            obj_node_v = obj_node_v->next;
            continue;
        }

        //-----------------------------

        // Weight checking.
        if (w_tag == 1)
            //            weight_new = weight > obj_v->weight ? weight : obj_v->weight;
            weight_new = weight;
        else //w_tag ==2
        {
            weight_new = weight + obj_v->weight;
            //			weight_new = comp_weight(w_tag, G) + obj_v->weight;

            if (foundResult && weight_new >= curWeight) {
                obj_node_v = obj_node_v->next;
                continue;
            }
        }
        //-----------------------------

        // Update the IF_S.
        bst_node_list_v = update_IF_obj(IF_S, obj_v);

        // Update the S_0.
        // obj_v is added at the first place of S_0.
        add_obj_set_entry(obj_v, G);

        // Recursive call
        tag = search_sub(m_opt, B, IF_S, G, q, d_o_q_new, weight_new, curBest, curWeight, foundResult);

        // Restore the S_0.
        remove_obj_set_entry(G);

        // Restore the IF_S.
        restore_IF_bst_node_list(IF_S, bst_node_list_v);
        release_bst_node_list(bst_node_list_v);

        if (w_tag == 1 && tag) {
            return true;
        }

        // Try the next object candidate.
        obj_node_v = obj_node_v->next;
    }
    return false;
}

obj_set_t* baselineSearch(int m_opt, B_KEY_TYPE B, query_t* q, psi_t* psi_v, obj_set_t* S, obj_t* o_t, obj_set_t*& G_curBest, W_TYPE& curWeight)
{
    obj_set_t* G;
    bst_t* IF_S;
    loc_t* loc_v;
    B_KEY_TYPE d_o_q;

    IF_S = const_IF_sorted(S, psi_v);

    G = alloc_obj_set();
    add_obj_set_entry(o_t, G);

    loc_v = get_obj_loc(o_t);
    d_o_q = calc_dist_loc(loc_v, q->loc_v);
    release_loc(loc_v);

    //check whether this IF_S can find a solution
    search_sub2(m_opt, B, IF_S, G, q, G_curBest, curWeight);

    release_IF(IF_S);
    release_obj_set(G);

    return G_curBest;
}

/*
 *	The sub-procedure of function "EnumerateSubset".
 *  curSet and curCost store the current best set and cost, respectively.
 *  x = d(o1,o2)
 *  cost = cost of the current set (updateable by cost_new when Sum, SumMax)
 */
bool search_sub2(int m_opt, B_KEY_TYPE B, bst_t* IF_S, obj_set_t* G, query_t* q, obj_set_t*& curBest, W_TYPE& curWeight)
{
    obj_t* obj_v;
    bst_node_t* bst_node_v;
    obj_node_t* obj_node_v;
    bst_node_list_t* bst_node_list_v;
    W_TYPE weight_new;
    bool tag;

    // if psi==\emptyset
    if (IF_S->node_n == 0) {
        //		printf("node_n==0\n");
        //update current best
        weight_new = comp_weight(w_tag, G);
        if (weight_new < curWeight) {
            release_obj_set(curBest);
            curBest = copy_obj_set(G);
            curWeight = weight_new;
        }
        return true;
    }

    bst_node_v = IF_S->root;
    obj_node_v = bst_node_v->p_list_obj->head->next;
    while (obj_node_v != NULL) { // Pick an object.
        obj_v = obj_node_v->obj_v;

        add_obj_set_entry(obj_v, G);

        // Distance checking.
        if (comp_cost(cost_tag, G, q) > B) {
            remove_obj_set_entry(G);
            obj_node_v = obj_node_v->next;
            continue;
        }

        // Weight checking.
        if (comp_weight(w_tag, G) >= curWeight) {
            remove_obj_set_entry(G);
            obj_node_v = obj_node_v->next;
            continue;
        }

        // Update the IF_S.
        bst_node_list_v = update_IF_obj(IF_S, obj_v);

        // Recursive call
        tag = search_sub2(m_opt, B, IF_S, G, q, curBest, curWeight);

        // Restore the S_0.
        remove_obj_set_entry(G);

        // Restore the IF_S.
        restore_IF_bst_node_list(IF_S, bst_node_list_v);
        release_bst_node_list(bst_node_list_v);

        //		if(w_tag==1 &&tag)
        //			return true;

        // Try the next object candidate.
        obj_node_v = obj_node_v->next;
    }
    return false;
}


obj_set_t* APPRO2(query_t* q, obj_set_t* S)
{
    psi_t* psi_v;
    obj_node_t *obj_node_v, *obj_node_picked;
    obj_set_t* obj_set_v;
    double ratio, ratio_max;
    k_node_t* k_node_v;
    int num;

    // Obtain the "un-covered" keywords by S.
    k_node_v = (k_node_t*)malloc(sizeof(k_node_t));
    memset(k_node_v, 0, sizeof(k_node_t));
    copy_k_list(k_node_v, q->psi_v->k_head);
    psi_v = const_psi(k_node_v);

    //    	printf("psi_v:%d\n",psi_v->key_n);
    //    	print_k_list(psi_v->k_head, stdout);
    //
    //	print_obj_set(S, stdout);

    //---prechecking
    int cnt = 0;
    obj_node_v = S->head->next;

    while (obj_node_v != NULL) {
        obj_node_v->numOfInt = is_relevant_obj(obj_node_v->obj_v, psi_v);
        cnt += obj_node_v->numOfInt;

        obj_node_v = obj_node_v->next;
    }
    //trivial case: each keyword is covered by one object only, return S
    if (cnt == psi_v->key_n) {
        release_psi(psi_v);
        obj_set_v = copy_obj_set(S);
        return obj_set_v;
    }
    //	//---

    obj_set_v = alloc_obj_set();
    // while there are keyword not covered yet
    while (psi_v->key_n > 0) {
        obj_node_picked = NULL;
        obj_node_v = S->head->next;

        ratio_max = 0.0;
        while (obj_node_v != NULL) {

            ratio = ((double)is_relevant_obj(obj_node_v->obj_v, psi_v)) / ((double)obj_node_v->obj_v->weight);
            if (ratio > ratio_max) {
                ratio_max = ratio;
                obj_node_picked = obj_node_v;
            }

            obj_node_v = obj_node_v->next;
        }
        if (obj_node_picked == NULL) {
            release_obj_set(obj_set_v);
            obj_set_v = NULL;
            break;
        }

        // update the uncovered keywords
        // S=S \cup o'
        // \psi = \psi - o'\psi
        add_obj_set_entry(obj_node_picked->obj_v, obj_set_v);
        psi_exclusion(psi_v, obj_node_picked->obj_v);
    }
    if (psi_v->key_n != 0) {
        release_obj_set(obj_set_v);
        obj_set_v = NULL;
    }
    // Release the memory.
    release_psi(psi_v);

    return obj_set_v;
}


/*
 * keep add objects with min weight until feasible
 * pre-checking: add all objects with weight <= o_w
 */
bool APPRO1_addver2(query_t* q, obj_set_t* S, obj_set_t* G, bst_t* Tab)
{

    obj_node_t* obj_node_v;
    bool feasible;

    b_heap_t* R_heap;
    int top;
    obj_t *o, *o_next;

    feasible = false;

    //	printf("S->obj_n:%d\n",S->obj_n);
    //Sort the objects in R by their weight
    R_heap = heap_sort_obj_set_weight(S);
    top = b_h_get_top(R_heap);

    o_next = R_heap->obj_arr[top].obj_v;

    while (true) {
        if (o_next == NULL)
            break;

        o = o_next;

        top = b_h_get_top(R_heap);
        if (top == 0)
            o_next = NULL;
        else
            o_next = R_heap->obj_arr[top].obj_v;

        // update Tab
        add_obj_set_entry(o, G);
        update_Tab(Tab->root, o);

        //print_bst(Tab);
        if (check_Tab_not_zero(Tab->root)) {
            feasible = true;
            break;
        }
    }

    release_b_heap(R_heap);
    //------------------------
    return feasible;
}


/*
 * Construct the Tab with keyword in @psi_v
 */
bst_t* build_Tab(psi_t* psi_v, obj_set_t* obj_set_v)
{
    bst_t* bst_v;
    bst_node_t* bst_node_v;
    obj_node_t* obj_node_v;

    bst_v = bst_ini();

    k_node_t* k_node_v = psi_v->k_head->next;
    while (k_node_v != NULL) {
        // Insert a new keyword in the binary tree.
        bst_node_v = (bst_node_t*)malloc(sizeof(bst_node_t));
        memset(bst_node_v, 0, sizeof(bst_node_t));

        /*s*/
        stat_v.memory_v += sizeof(bst_node_t);
        if (stat_v.memory_v > stat_v.memory_max)
            stat_v.memory_max = stat_v.memory_v;
        /*s*/

        // Update the posting list.
        bst_node_v->key = k_node_v->key;
        bst_node_v->freq = 0;

        //---
        // count the frequency of the keyword
        obj_node_v = obj_set_v->head->next;
        while (obj_node_v != NULL) {

            if (has_key_obj(obj_node_v->obj_v, k_node_v->key))
                bst_node_v->freq++;

            obj_node_v = obj_node_v->next;
        }
        //---

        bst_insert(bst_v, bst_node_v);
        k_node_v = k_node_v->next;
    }

    return bst_v;
}


/*
 * Construct the Tab with keyword in @psi_v with min weight|<=seed_weight object only
 */
bst_t* build_Tab_min_weight(psi_t* psi_v, obj_set_t* obj_set_v, W_TYPE seed_weight)
{
	bst_t* bst_v;
	bst_node_t* bst_node_v;
	obj_node_t* obj_node_v;

	bst_v = bst_ini();

	k_node_t* k_node_v = psi_v->k_head->next;
	while (k_node_v != NULL) {
		// Insert a new keyword in the binary tree.
		bst_node_v = (bst_node_t*)malloc(sizeof(bst_node_t));
		memset(bst_node_v, 0, sizeof(bst_node_t));

		/*s*/
		stat_v.memory_v += sizeof(bst_node_t);
		if (stat_v.memory_v > stat_v.memory_max)
			stat_v.memory_max = stat_v.memory_v;
		/*s*/

		// Update the posting list.
		bst_node_v->key = k_node_v->key;
		bst_node_v->freq = 0;

		//---
		obj_t* obj_v1 = NULL; //min weight one
		W_TYPE min_weight = MAX_WEIGHT;

		// count the frequency of the keyword
		obj_node_v = obj_set_v->head->next;
		while (obj_node_v != NULL) {

			if (has_key_obj(obj_node_v->obj_v, k_node_v->key) && obj_node_v->obj_v->weight < min_weight){
				
//				bst_node_v->freq++;
				min_weight = obj_node_v->obj_v->weight;
				obj_v1 = obj_node_v->obj_v;
					if(min_weight<=seed_weight)
						break;
			}
			obj_node_v = obj_node_v->next;
		}
		
		bst_node_v->obj_v1 = obj_v1;
		bst_node_v->min_weight = min_weight;
		//---

		bst_insert(bst_v, bst_node_v);
		k_node_v = k_node_v->next;
	}

	return bst_v;
}


// return true if all node freq >0
bool get_min_weight_set(bst_node_t* bst_node_v, obj_set_t* obj_set_v)
{
	if (bst_node_v == NULL)
		return true;
	
	if(bst_node_v->obj_v1==NULL)
		return false;
	
	//cur is ok
	add_obj_set_entry(bst_node_v->obj_v1, obj_set_v);
		
	return get_min_weight_set(bst_node_v->left, obj_set_v) && get_min_weight_set(bst_node_v->right, obj_set_v);
}

// return true if all node freq >0
bool check_Tab_not_zero(bst_node_t* bst_node_v)
{
    if (bst_node_v == NULL)
        return true;
    return bst_node_v->freq > 0 && check_Tab_not_zero(bst_node_v->left) && check_Tab_not_zero(bst_node_v->right);
}

// return true if all node freq ==1
bool check_Tab_all_one(bst_node_t* bst_node_v)
{
    if (bst_node_v == NULL)
        return true;
    return bst_node_v->freq == 1 && check_Tab_all_one(bst_node_v->left) && check_Tab_all_one(bst_node_v->right);
}

//return true if obj_v is a critical obj - (at least) one of the freq will become 0
bool check_critical_obj_in_Tab(bst_node_t* bst_node_v, obj_t* obj_v)
{
    if (bst_node_v == NULL)
        return false;
    bool tag = false;

    if (has_key_obj(obj_v, bst_node_v->key))
        bst_node_v->freq--;

    if (bst_node_v->freq == 0)
        return true;

    //cannot use || here because it may skip the second part (the freq--)
    if (check_critical_obj_in_Tab(bst_node_v->left, obj_v))
        tag = true;

    if (check_critical_obj_in_Tab(bst_node_v->right, obj_v))
        tag = true;

    return tag;
}

void update_Tab(bst_node_t* bst_node_v, obj_t* obj_v)
{
    if (bst_node_v == NULL)
        return;

    if (has_key_obj(obj_v, bst_node_v->key))
        bst_node_v->freq++;

    update_Tab(bst_node_v->left, obj_v);
    update_Tab(bst_node_v->right, obj_v);

    return;
}

/*
 *	Sort the objects in @obj_set_v by their weight.
 *
 *	Method: heap-sort.
 *	Alternative method: based on the binary search tree.
 */
b_heap_t* heap_sort_obj_set_weight(obj_set_t* obj_set_v)
{
    int cur;
    B_KEY_TYPE dist;
    b_heap_t* b_h;
    obj_node_t* obj_node_v;
    loc_t* loc_v;

    b_h = alloc_b_heap(obj_set_v->obj_n + 1);

    cur = 1;
    obj_node_v = obj_set_v->head->next;
    while (obj_node_v != NULL) {
        b_h->obj_arr[cur].key = obj_node_v->obj_v->weight;
        b_h->obj_arr[cur].obj_v = obj_node_v->obj_v;

        b_h_insert(b_h, cur);

        cur++;

        obj_node_v = obj_node_v->next;
    }

    return b_h;
}


/*
 * Find the most infrequent query keyword in inverted list index @IF_bst_v
 *
 * Return the node with the most infrequent query keyword
 *
 * return NULL if no keyword is found in @IF_bst_v
 *
 */
bst_node_t* findMostInfreqKey(bst_t* IF_bst_v, psi_t* psi_v)
{
    bst_node_t* mostInfreqKey = NULL;
    B_KEY_TYPE count = INFINITY;
    k_node_t* k_node_v = psi_v->k_head->next;
    bst_node_t* bst_node_v;

    while (k_node_v != NULL) {
        bst_node_v = bst_search(IF_bst_v, k_node_v->key);

        if (bst_node_v != NULL && count > bst_node_v->p_list_obj->obj_n) {
            mostInfreqKey = bst_node_v;
            count = bst_node_v->p_list_obj->obj_n;
        }

        k_node_v = k_node_v->next;
    }

    // printf("key:%f \t cnt:%f\n",mostInfreqKey, count);

    return mostInfreqKey;
}

W_TYPE comp_weight(int weight_tag, obj_set_t* obj_set_v)
{
    if (weight_tag == 1)
        return comp_weight_max(obj_set_v);
    else
        return comp_weight_sum(obj_set_v);
}

W_TYPE comp_weight_max(obj_set_t* obj_set_v)
{
    obj_node_t* obj_node_v;
    W_TYPE maxWeight = 0;

    obj_node_v = obj_set_v->head->next;
    while (obj_node_v != NULL) {
        if (obj_node_v->obj_v->weight > maxWeight)
            maxWeight = obj_node_v->obj_v->weight;
        obj_node_v = obj_node_v->next;
    }
    return maxWeight;
}

W_TYPE comp_weight_sum(obj_set_t* obj_set_v)
{
    obj_node_t* obj_node_v;
    W_TYPE sumWeight = 0;

    obj_node_v = obj_set_v->head->next;
    while (obj_node_v != NULL) {
        sumWeight += obj_node_v->obj_v->weight;
        obj_node_v = obj_node_v->next;
    }
    return sumWeight;
}

//--------------------------------------------------------------------------------------
// IR-tree range query related

/*
 *	Retrieve all the objects located at the sub-tree rooted at @node_v.
 *	The retrieved objects are stored in obj_set_v.
 */

void retrieve_sub_tree_within_weight(node_t* node_v, obj_set_t*& obj_set_v, query_t* q,
    W_TYPE max_weight)
{
    int i;
    BIT_TYPE p_list;

    if (node_v->level == 0) {
        // node_v is a leaf-node.
        // Retrieve all its objects.
        for (i = 0; i < node_v->num; i++) {
            if (is_relevant_obj((obj_t*)(node_v->child[i]), q)
                && ((obj_t*)(node_v->child[i]))->weight <= max_weight)
                add_obj_set_entry((obj_t*)(node_v->child[i]), obj_set_v);
            //                add_obj_set_entry_sorted((obj_t*)(node_v->child[i]), obj_set_v, ((obj_t*)(node_v->child[i]))->weight);
        }
    } else {
        // node_v is an inner-node.
        // Invoke the function recursively.
        p_list = is_relevant_node_within_weight(node_v, q, max_weight);
        for (i = 0; i < node_v->num; i++) {
            if (get_k_bit(p_list, i))
                retrieve_sub_tree_within_weight((node_t*)(node_v->child[i]), obj_set_v, q,
                    max_weight);
        }
    }
}

//===========================================================================================
// range query for intersection of two circle (disk) region

/*
 *	Range query on the sub-tree rooted at @node_v.
 *	@disk_v indicates the range which is a circle.
 *
 *	The results are stored in @obj_set_v.
 */
void range_query_sub(node_t* node_v, disk_t* disk_v1, disk_t* disk_v2, obj_set_t*& obj_set_v, query_t* q)
{
	int i;
	BIT_TYPE p_list;
	range* MBR;

	if (node_v->parent == NULL)
		MBR = get_MBR_node(node_v, IRTree_v.dim);
	else
		MBR = node_v->parent->MBRs[node_v->loc];

	//No overlapping.
	if (!is_overlap(MBR, disk_v1) || !is_overlap(MBR, disk_v2))
		return;

	//Enclosed entrely.
	if (is_enclosed(MBR, disk_v1) && is_enclosed(MBR, disk_v2)) {
		retrieve_sub_tree(node_v, obj_set_v, q);
		if (node_v->parent == NULL) {
			free(MBR);

			/*s*/
			stat_v.memory_v -= IRTree_v.dim * sizeof(range);
			/*s*/
		}

		return;
	}

	//The remaining cases.
	if (node_v->level == 0) //node_v is a leaf-node.
	{
		///for each object inside the leaf node
		for (i = 0; i < node_v->num; i++) {
			///if inside the disk and relevant, then add into obj_set_v
			if (is_enclosed(((obj_t*)(node_v->child[i]))->MBR, disk_v1) && is_enclosed(((obj_t*)(node_v->child[i]))->MBR, disk_v2) && is_relevant_obj((obj_t*)(node_v->child[i]), q))
				add_obj_set_entry((obj_t*)(node_v->child[i]), obj_set_v);
		}
	} else //node_v is an inner-node.
	{
		///retrieve a list of relevant childean
		p_list = is_relevant_node(node_v, q);

		///recursive call for each child in the list
		for (i = 0; i < node_v->num; i++) {
			if (get_k_bit(p_list, i))
				range_query_sub((node_t*)(node_v->child[i]), disk_v1, disk_v2, obj_set_v, q);
		}
	}

	if (node_v->parent == NULL) {
		free(MBR);
		/*s*/
		stat_v.memory_v -= IRTree_v.dim * sizeof(range);
		/*s*/
	}
}

/*
 *	Circle range query.
 *
 *	DFS: recursive implementation.
 */
obj_set_t* range_query(disk_t* disk_v1, disk_t* disk_v2, query_t* q)
{
	obj_set_t* obj_set_v;

	obj_set_v = alloc_obj_set();
	range_query_sub(IRTree_v.root, disk_v1, disk_v2, obj_set_v, q);

	return obj_set_v;
}
//===========================================================================================
/*
 *	Range query on the sub-tree rooted at @node_v.
 *	@disk_v indicates the range which is a circle.
 *
 *	The results are stored in @obj_set_v.
 */
void range_query_within_weight_sub(node_t* node_v, disk_t* disk_v1, disk_t* disk_v2, obj_set_t*& obj_set_v, query_t* q, W_TYPE max_weight)
{
    int i;
    BIT_TYPE p_list;
    range* MBR;

    if (node_v->parent == NULL)
        MBR = get_MBR_node(node_v, IRTree_v.dim);
    else
        MBR = node_v->parent->MBRs[node_v->loc];

    // No overlapping.
    if (!is_overlap(MBR, disk_v1))
        return;

    if (disk_v2 != NULL && !is_overlap(MBR, disk_v2))
        return;

    // Enclosed entrely.
    if (is_enclosed(MBR, disk_v1) && (disk_v2 == NULL || is_enclosed(MBR, disk_v2))) {
        retrieve_sub_tree_within_weight(node_v, obj_set_v, q, max_weight);
        if (node_v->parent == NULL) {
            free(MBR);

            /*s*/
            stat_v.memory_v -= IRTree_v.dim * sizeof(range);
            /*s*/
        }

        return;
    }

    // The remaining cases.
    if (node_v->level == 0) // node_v is a leaf-node.
    {
        /// for each object inside the leaf node
        for (i = 0; i < node_v->num; i++) {
            /// if inside the disk and relevant, then add into obj_set_v
            if (is_enclosed(((obj_t*)(node_v->child[i]))->MBR, disk_v1)
                && (disk_v2 == NULL || is_enclosed(((obj_t*)(node_v->child[i]))->MBR, disk_v2))
                && is_relevant_obj((obj_t*)(node_v->child[i]), q)
                && ((obj_t*)(node_v->child[i]))->weight <= max_weight)
                add_obj_set_entry((obj_t*)(node_v->child[i]), obj_set_v);
            //                add_obj_set_entry_sorted((obj_t*)(node_v->child[i]), obj_set_v, ((obj_t*)(node_v->child[i]))->weight);
        }
    } else // node_v is an inner-node.
    {
        /// retrieve a list of relevant childean
        p_list = is_relevant_node_within_weight(node_v, q, max_weight);

        /// recursive call for each child in the list
        for (i = 0; i < node_v->num; i++) {
            if (get_k_bit(p_list, i))
                range_query_within_weight_sub((node_t*)(node_v->child[i]), disk_v1, disk_v2,
                    obj_set_v, q, max_weight);
        }
    }

    if (node_v->parent == NULL) {
        free(MBR);
        /*s*/
        stat_v.memory_v -= IRTree_v.dim * sizeof(range);
        /*s*/
    }
}

/*
 *	Circle range query.
 *
 *	DFS: recursive implementation.
 *
 * only retreive objects with < @max_weight
 * **not <=
 */
obj_set_t* range_query_within_weight(disk_t* disk_v1, disk_t* disk_v2, query_t* q, W_TYPE max_weight)
{
    obj_set_t* obj_set_v;

    obj_set_v = alloc_obj_set();
    range_query_within_weight_sub(IRTree_v.root, disk_v1, disk_v2, obj_set_v, q, max_weight);

    return obj_set_v;
}

//--------------------------------------------------------------------------------------
// IR-tree range query related 2

/*
 *	Retrieve all the objects located at the sub-tree rooted at @node_v.
 *	The retrieved objects are stored in obj_set_v.
 */

void retrieve_sub_tree_ellipse_within_weight(node_t* node_v, obj_set_t*& obj_set_v, query_t* q,
    W_TYPE max_weight)
{
    int i;
    BIT_TYPE p_list;

    if (node_v->level == 0) {
        // node_v is a leaf-node.
        // Retrieve all its objects.
        for (i = 0; i < node_v->num; i++) {
            if (is_relevant_obj((obj_t*)(node_v->child[i]), q)
                && ((obj_t*)(node_v->child[i]))->weight < max_weight)
                add_obj_set_entry((obj_t*)(node_v->child[i]), obj_set_v);
            //                add_obj_set_entry_sorted((obj_t*)(node_v->child[i]), obj_set_v, ((obj_t*)(node_v->child[i]))->weight);
        }
    } else {
        // node_v is an inner-node.
        // Invoke the function recursively.
        p_list = is_relevant_node_within_weight(node_v, q, max_weight);
        for (i = 0; i < node_v->num; i++) {
            if (get_k_bit(p_list, i))
                retrieve_sub_tree_ellipse_within_weight((node_t*)(node_v->child[i]), obj_set_v, q,
                    max_weight);
        }
    }
}

/*
 *	Range query on the sub-tree rooted at @node_v.
 *	@disk_v indicates the range which is a circle.
 *
 *	The results are stored in @obj_set_v.
 */
void range_query_ellipse_within_weight_sub(node_t* node_v, disk_t* disk_v1, obj_set_t*& obj_set_v, query_t* q, W_TYPE max_weight, ellipse_t* ellipse_v)
{
    int i;
    BIT_TYPE p_list;
    range* MBR;

    if (node_v->parent == NULL)
        MBR = get_MBR_node(node_v, IRTree_v.dim);
    else
        MBR = node_v->parent->MBRs[node_v->loc];

    // No overlapping.
    if (!is_overlap(MBR, disk_v1))
        return;

    //---
    //    if (!is_overlap(MBR, ellipse_v))
    //		return;
    //---

    // Enclosed entrely.
    if (is_enclosed(MBR, disk_v1) && is_enclosed(MBR, ellipse_v)) {
        retrieve_sub_tree_ellipse_within_weight(node_v, obj_set_v, q, max_weight);
        if (node_v->parent == NULL) {
            free(MBR);

            /*s*/
            stat_v.memory_v -= IRTree_v.dim * sizeof(range);
            /*s*/
        }

        return;
    }

    // The remaining cases.
    if (node_v->level == 0) // node_v is a leaf-node.
    {
        /// for each object inside the leaf node
        for (i = 0; i < node_v->num; i++) {
            /// if inside the disk and relevant, then add into obj_set_v
            if (is_enclosed(((obj_t*)(node_v->child[i]))->MBR, disk_v1)
                && is_relevant_obj((obj_t*)(node_v->child[i]), q)
                && ((obj_t*)(node_v->child[i]))->weight < max_weight
                && is_enclosed(((obj_t*)(node_v->child[i]))->MBR, ellipse_v))
                add_obj_set_entry((obj_t*)(node_v->child[i]), obj_set_v);
            //                add_obj_set_entry_sorted((obj_t*)(node_v->child[i]), obj_set_v, ((obj_t*)(node_v->child[i]))->weight);
        }
    } else // node_v is an inner-node.
    {
        /// retrieve a list of relevant childean
        p_list = is_relevant_node_within_weight(node_v, q, max_weight);

        /// recursive call for each child in the list
        for (i = 0; i < node_v->num; i++) {
            if (get_k_bit(p_list, i))
                range_query_ellipse_within_weight_sub((node_t*)(node_v->child[i]), disk_v1,
                    obj_set_v, q, max_weight, ellipse_v);
        }
    }

    if (node_v->parent == NULL) {
        free(MBR);
        /*s*/
        stat_v.memory_v -= IRTree_v.dim * sizeof(range);
        /*s*/
    }
}

/*
 *	Circle range query.
 *
 *	DFS: recursive implementation.
 *
 * only retreive objects with < @max_weight
 *
 * utilize the ellipse pruning
 */
obj_set_t* range_query_ellipse_within_weight(disk_t* disk_v1, query_t* q, W_TYPE max_weight, ellipse_t* ellipse_v)
{
    obj_set_t* obj_set_v;

    obj_set_v = alloc_obj_set();
    range_query_ellipse_within_weight_sub(IRTree_v.root, disk_v1, obj_set_v, q, max_weight, ellipse_v);

    return obj_set_v;
}

//-----------------------------------------------------------------------
/*
 *	Check whether the sub-tree rooted at a node @node_v
 *	contains a keyword @key.
 *
 */
BIT_TYPE has_key_node_within_weight(node_t* node_v, KEY_TYPE key, W_TYPE max_weight)
{
    bst_node_t* bst_node_v;

    if ((bst_node_v = bst_search(node_v->bst_v, key))) {
        // return bst_node_v->p_list;
        //----------------
        // perform additional checking/filtering on weight
        // if weight of the node is > max_weight
        // set the bit to 0, i.e., del k bit
        BIT_TYPE p_list = bst_node_v->p_list;
        for (int i = 0; i < node_v->num; i++) {
            if (get_k_bit(p_list, i) && bst_node_v->p_list_weight[i] > max_weight)
                delete_k_bit(p_list, i);
        }
        //----------------
        return p_list;
    }
    return 0;
}

/*
 *	Check whether a node @node_v is "relevant" or not.
 *
 **  check whether the children of this @node_v contain any keyword
 ** return all children's positions that contain query keyword by 32-bits (e.g.,
 *1001)
 */
BIT_TYPE is_relevant_node_within_weight(node_t* node_v, query_t* q, W_TYPE weight)
{
    BIT_TYPE res, res_t;

    k_node_t* k_node_v;

    res = 0;
    k_node_v = q->psi_v->k_head->next;
    while (k_node_v != NULL) {
        res_t = has_key_node_within_weight(node_v, k_node_v->key, weight);
        union_bit(res, res_t);

        if (res == UINT_MAX)
            return res;

        k_node_v = k_node_v->next;
    }

    return res;
}

//--------------------------------------------------------------------------------------

/*
 *	Construct the IF on a set of objects @obj_set_v for the keywords in
 *@psi_v.
 *
 *	1. The IF structure is indexed by a binary search tree.
 *	2. No ordering is imposed in IF.
 *  3. Each list is ordered by object weight in ascending order
 */
bst_t* const_IF_sorted(obj_set_t* obj_set_v, psi_t* psi_v)
{
    int i;
    bst_t* IF_bst_v;
    k_node_t* k_node_v;
    obj_node_t* obj_node_v;
    bst_node_t* bst_node_v;

    IF_bst_v = bst_ini();

    k_node_v = psi_v->k_head->next;
    /// for each keyword in psi_v
    for (i = 0; i < psi_v->key_n; i++) {

        bst_node_v = (bst_node_t*)malloc(sizeof(bst_node_t));
        memset(bst_node_v, 0, sizeof(bst_node_t));

        /*s*/
        stat_v.memory_v += sizeof(bst_node_t);
        /*s*/

        bst_node_v->key = k_node_v->key;
        bst_node_v->p_list_obj = alloc_obj_set();

        /// for each objects in obj_set_v
        obj_node_v = obj_set_v->head->next;
        while (obj_node_v != NULL) {
            if (has_key_obj(obj_node_v->obj_v, k_node_v->key)) {
                add_obj_set_entry_sorted(obj_node_v->obj_v, bst_node_v->p_list_obj,
                    obj_node_v->obj_v->weight);
            }
            obj_node_v = obj_node_v->next;
        }

        bst_insert(IF_bst_v, bst_node_v);

        k_node_v = k_node_v->next;
    }

    /*s*/
    if (stat_v.memory_v > stat_v.memory_max)
        stat_v.memory_max = stat_v.memory_v;
    /*s*/

    return IF_bst_v;
}

/*
 * Insert the object @obj_v into the inverted file @IF_bst_v
 * Note that only keywords in @IF_bst_v are considered
 *
 */
void insert_IF_obj(bst_t* IF_bst_v, obj_t* obj_v)
{
    k_node_t* k_node_v;
    bst_node_t* bst_node_v;

    k_node_v = obj_v->k_head->next;
    while (k_node_v != NULL) {
        if ((bst_node_v = bst_search(IF_bst_v, k_node_v->key)) != NULL) {
            add_obj_set_entry(obj_v, bst_node_v->p_list_obj);
        }
        k_node_v = k_node_v->next;
    }

    /*s*/
    if (stat_v.memory_v > stat_v.memory_max)
        stat_v.memory_max = stat_v.memory_v;
    /*s*/

    return;
}

// Function to find N-th Harmonic Number
double nthHarmonic(int N)
{
    // H1 = 1
    float harmonic = 1.00;

    // loop to apply the forumula
    // Hn = H1 + H2 + H3 ... + Hn-1 + Hn-1 + 1/n
    for (int i = 2; i <= N; i++) {
        harmonic += (float)1 / i;
    }

    return harmonic;
}

/*
 *	Return true if the IF @bst_v contain at least 1 obj in each inverted list
 */
bool IF_obj_check(bst_node_t* bst_node_v)
{
    if (bst_node_v == NULL)
        return true;

    if (bst_node_v->p_list_obj == NULL || bst_node_v->p_list_obj->obj_n == 0)
        return false;

    if (IF_obj_check(bst_node_v->left) && IF_obj_check(bst_node_v->right))
        return true;

    return false;
}

/*
 *	For each query keyword, we find the object with minimum weight
 *	The maximum among these objects correspond to the minimum weight to cover all keywords
 *
 */
W_TYPE findOptMaxWeightLB(bst_t* IF_bst_v, query_t* q, B_KEY_TYPE B)
{
    W_TYPE maxWeight;
    W_TYPE minAmongObj;

    k_node_t* k_node_v;
    bst_node_t* bst_node_v;
    obj_node_t* obj_node_v;

    maxWeight = -1;

    k_node_v = q->psi_v->k_head->next;
    while (k_node_v != NULL) {
        minAmongObj = MAX_WEIGHT;

        bst_node_v = bst_search(IF_bst_v, k_node_v->key);
        //		printf("key:%.0lf\n",k_node_v->key);
        //		print_obj_set(bst_node_v->p_list_obj, stdout);

        obj_node_v = bst_node_v->p_list_obj->head->next;
        while (obj_node_v != NULL) {
            if (obj_node_v->obj_v->weight < minAmongObj && calc_minDist(obj_node_v->obj_v->MBR, q->loc_v) <= B)
                minAmongObj = obj_node_v->obj_v->weight;
            obj_node_v = obj_node_v->next;
        }
        if (minAmongObj > maxWeight)
            maxWeight = minAmongObj;
        //		printf("maxWeight:%d\n",maxWeight);
        k_node_v = k_node_v->next;
    }
    return maxWeight;
}

/*
 *	For each query keyword, we find the object with minimum weight
 *	The set of these objects form a set
 *
 * 	for the object with same (min) weight
 * 	the one with min d(o,q) is picked
 *
 */
obj_set_t* findMinWeightSet(bst_t* IF_bst_v, query_t* q, B_KEY_TYPE B)
{
    obj_set_t* obj_set_v;
    W_TYPE minAmongObj;

    B_KEY_TYPE dist, minAmongObjDist;

    k_node_t* k_node_v;
    bst_node_t* bst_node_v;
    obj_node_t *obj_node_v, *obj_node_min;

    obj_set_v = alloc_obj_set();

    k_node_v = q->psi_v->k_head->next;
    while (k_node_v != NULL) {
        minAmongObj = MAX_WEIGHT;
        minAmongObjDist = B;
        obj_node_min = NULL;
        bst_node_v = bst_search(IF_bst_v, k_node_v->key);
        //		printf("key:%.0lf\n",k_node_v->key);
        //		print_obj_set(bst_node_v->p_list_obj, stdout);

        obj_node_v = bst_node_v->p_list_obj->head->next;
        while (obj_node_v != NULL) {
            dist = calc_minDist(obj_node_v->obj_v->MBR, q->loc_v);
            if (dist <= B) {
                if (obj_node_v->obj_v->weight == minAmongObj) {
                    if (dist < minAmongObjDist) {
                        minAmongObj = obj_node_v->obj_v->weight;
                        minAmongObjDist = dist;
                        obj_node_min = obj_node_v;
                    }
                } else if (obj_node_v->obj_v->weight < minAmongObj) {
                    minAmongObj = obj_node_v->obj_v->weight;
                    minAmongObjDist = dist;
                    obj_node_min = obj_node_v;
                }
            }
            obj_node_v = obj_node_v->next;
        }

        //---
        if (obj_node_min != NULL)
            add_obj_set_entry(obj_node_min->obj_v, obj_set_v);
        else {
            //do not exist an object in D(q,B) cover this key
            release_obj_set(obj_set_v);
            return NULL;
        }
        //---
        //		printf("maxWeight:%d\n",maxWeight);
        k_node_v = k_node_v->next;
    }
    return obj_set_v;
}

/*
 *	For each query keyword, we find the object with minimum "average" weight
 *	The sum among these objects correspond to the minimum weight to cover all keywords
 *
 */
W_TYPE findOptSumWeightLB(bst_t* IF_bst_v, query_t* q, B_KEY_TYPE B)
{
    B_KEY_TYPE sumWeight;
    B_KEY_TYPE minAvgWeight;
    B_KEY_TYPE dist, tmp;

    B_KEY_TYPE sumNumOfObj;
    B_KEY_TYPE minNumOfObj;

    k_node_t* k_node_v;
    bst_node_t* bst_node_v;
    obj_node_t* obj_node_v;

    sumWeight = 0.0;
    sumNumOfObj = 0.0;

    k_node_v = q->psi_v->k_head->next;
    while (k_node_v != NULL) {
        minAvgWeight = MAX_WEIGHT;
        minNumOfObj = 2;

        bst_node_v = bst_search(IF_bst_v, k_node_v->key);
        //		printf("key:%.0lf\n",k_node_v->key);
        //		print_obj_set(bst_node_v->p_list_obj, stdout);

        obj_node_v = bst_node_v->p_list_obj->head->next;
        while (obj_node_v != NULL) {
            dist = calc_minDist(obj_node_v->obj_v->MBR, q->loc_v);
            if (dist <= B) {

                int numInter = number_intersection(obj_node_v->obj_v->k_head, q->psi_v->k_head);
                if (numInter > 0) {
                    tmp = obj_node_v->obj_v->weight / numInter;
                    if (tmp < minAvgWeight)
                        minAvgWeight = tmp;
                    if (1 / numInter < minNumOfObj)
                        minNumOfObj = 1 / numInter;
                }
            }
            obj_node_v = obj_node_v->next;
        }

        if (minNumOfObj == 2) {
            ///not feasible within B
            return INFINITY;
        }

        sumWeight += minAvgWeight;
        sumNumOfObj += minNumOfObj;
        //		printf("maxWeight:%d\n",maxWeight);
        k_node_v = k_node_v->next;
    }
    //	printf("sumWeight:%lf\n",sumWeight);
    //	printf("sumNumOfObj:%lf\n",sumNumOfObj);

    return sumWeight;
}

void print_IF(bst_t* T)
{
    print_IF_node(T->root);
    printf("\n");
}

void print_IF_node(bst_node_t* bst_node_v)
{
    if (bst_node_v != NULL) {

        printf("key:%lf\t", bst_node_v->key);
        if (bst_node_v->p_list_obj != NULL) {
            obj_node_t* obj_node_v = bst_node_v->p_list_obj->head->next;
            while (obj_node_v != NULL) {
                printf("%d ", obj_node_v->obj_v->id);
                obj_node_v = obj_node_v->next;
            }
            printf("\n");
        }
        print_IF_node(bst_node_v->left);
        print_IF_node(bst_node_v->right);
    }
}


obj_set_t* comp_LB(query_t* q, disk_t* disk_q_B, B_KEY_TYPE& distLB)
{
    B_KEY_TYPE LB, dist;
    k_node_t* k_node_v;

    obj_set_t* obj_set_v;
    obj_t* obj_v;
    loc_t* loc_v;

    obj_set_v = alloc_obj_set();

    //Compute LB.
    LB = 0;
    k_node_v = q->psi_v->k_head->next;
    while (k_node_v != NULL) {
        obj_v = const_NN_key(q->loc_v, k_node_v->key, disk_q_B);
        if (obj_v == NULL) {
            distLB = INFINITY;
            release_obj_set(obj_set_v);
            return NULL;
        }

        loc_v = get_obj_loc(obj_v);
        dist = calc_dist_loc(loc_v, q->loc_v);
        if (dist > LB)
            LB = dist;

        release_loc(loc_v);
        add_obj_set_entry(obj_v, obj_set_v);
        k_node_v = k_node_v->next;
    }

    distLB = LB;
    return obj_set_v;
}

/* check objects one by one: remove objects that are being dominated by other object
 if @returnRemoved = true, return those removed objects
 */
obj_set_t* removeObjDominated_new(obj_set_t* O_t, psi_t* psi_v, bool returnRemoved)
{
    obj_set_t* obj_set_v;
    obj_node_t *obj_node_v, *obj_node_v2, *obj_node_temp;
    bool remove;

    if (returnRemoved)
        obj_set_v = alloc_obj_set();

    obj_node_v = O_t->head;
    // checking obj_node_v->next
    while (obj_node_v->next != NULL) {
        remove = false;
        obj_node_v2 = O_t->head->next;
        while (obj_node_v2 != NULL) {
            //skip the same object
            if (obj_node_v->next->obj_v == obj_node_v2->obj_v) {
                obj_node_v2 = obj_node_v2->next;
                continue;
            }

            if (obj_node_v->next->obj_v->weight >= obj_node_v2->obj_v->weight
                && isDominated(obj_node_v->next->obj_v, obj_node_v2->obj_v, psi_v)) {
                //remove obj_node_v->next, i.e. outer loop obj
                obj_node_temp = obj_node_v->next;
                obj_node_v->next = obj_node_v->next->next;
                //---
                if (returnRemoved) {
                    obj_node_temp->next = obj_set_v->head->next;
                    obj_set_v->head->next = obj_node_temp;
                    obj_set_v->obj_n++;
                    //---
                } else {
                    free(obj_node_temp);
                    /*s*/
                    stat_v.memory_v -= sizeof(obj_node_t);
                    /*s*/
                }
                remove = true;
                O_t->obj_n--;
                break;
            }
            obj_node_v2 = obj_node_v2->next;
        }
        if (!remove)
            obj_node_v = obj_node_v->next;
    }

    if (returnRemoved)
        return obj_set_v;
    else
        return NULL;
}

/* faster implementaiton of finding dominate objects
 * maintain a set of dominate objects
 *
 */
obj_set_t* pickDominateObj(obj_set_t* O_t, psi_t* psi_v)
{
    obj_set_t* obj_set_v;
    obj_node_t *obj_node_v, *obj_node_v2, *obj_node_temp;
    bool dominatingObj;

    obj_set_v = alloc_obj_set();

    obj_node_v = O_t->head->next;
    while (obj_node_v != NULL) {
        dominatingObj = true;
        obj_node_v2 = obj_set_v->head->next;
        while (obj_node_v2 != NULL) {
            //check if the current object is dominated by v2 in @obj_set_v
            if (obj_node_v->obj_v->weight >= obj_node_v2->obj_v->weight
                && isDominated(obj_node_v->obj_v, obj_node_v2->obj_v, psi_v)) {
                dominatingObj = false;
                break;
            }
            obj_node_v2 = obj_node_v2->next;
        }

        if (dominatingObj) {
            //check if the current object dominate v2->next in @obj_set_v
            obj_node_v2 = obj_set_v->head;
            while (obj_node_v2->next != NULL) {
                if (obj_node_v2->next->obj_v->weight >= obj_node_v->obj_v->weight
                    && isDominated(obj_node_v2->next->obj_v, obj_node_v->obj_v, psi_v)) {
                    //--
                    //remove v2->next
                    obj_node_temp = obj_node_v2->next;
                    obj_node_v2->next = obj_node_v2->next->next;
                    free(obj_node_temp);
                    obj_set_v->obj_n--;
                    //--
                } else
                    obj_node_v2 = obj_node_v2->next;
            }
            //add current object to obj_set_v
            add_obj_set_entry(obj_node_v->obj_v, obj_set_v);
        }
        obj_node_v = obj_node_v->next;
    }

    return obj_set_v;
}

/*
 1. filter the objects with weight > max_weight
 i.e., removing from object set @S
 2. for each removed object, we put it in a bucket
 return the an array of object set, corresponds to buckets
 */
obj_set_t** filterAndDistributeObj(obj_set_t* O_t, psi_t* psi_v, W_TYPE max_weight)
{
    obj_set_t** obj_set_set;
    obj_node_t *obj_node_v, *obj_node_v2, *obj_node_temp;
    bool remove;

    //initialize the array
    obj_set_set = (obj_set_t**)malloc(sizeof(obj_set_t*) * 10);
    memset(obj_set_set, 0, sizeof(obj_set_t*) * 10);
    for (int i = 0; i < 10; i++)
        obj_set_set[i] = alloc_obj_set();

    obj_node_v = O_t->head;
    // checking obj_node_v->next
    while (obj_node_v->next != NULL) {
        remove = false;
        obj_node_v2 = O_t->head->next;

        if (obj_node_v->next->obj_v->weight > max_weight) {
            //remove obj_node_v->next,
            obj_node_temp = obj_node_v->next;
            obj_node_v->next = obj_node_v->next->next;
            //---
            //put it into the corresponding bucket
            W_TYPE bucketNum = ceil(obj_node_temp->obj_v->weight * 10 / (MAX_WEIGHT)) - 1;
            obj_node_temp->next = obj_set_set[bucketNum]->head->next;
            obj_set_set[bucketNum]->head->next = obj_node_temp;
            obj_set_set[bucketNum]->obj_n++;
            //---
            //				free(obj_node_temp);
            remove = true;
            O_t->obj_n--;
        }
        if (!remove)
            obj_node_v = obj_node_v->next;
    }

    return obj_set_set;
}

/*
 *	Check whether a MBR @MBR is enclosed entirely by an ellipse @disk_v.
 */
bool is_enclosed(range* MBR, ellipse_t* ellipse_v)
{

    if (calc_maxDist(MBR, ellipse_v->center1, ellipse_v->center2) <= ellipse_v->diam)
        return true;
    return false;
}

/*
 *	Calculate the maxDist between a MBR @MBR and
 */
B_KEY_TYPE calc_maxDist(range* MBR, loc_t* loc_v1, loc_t* loc_v2)
{
    B_KEY_TYPE maxDist, dist, dist1, dist2;

    maxDist = 0;

    dist = 0;

    dist1 = pow(MBR[0].max - loc_v1->coord[0], 2) + pow(MBR[1].max - loc_v1->coord[1], 2);
    dist2 = pow(MBR[0].max - loc_v2->coord[0], 2) + pow(MBR[1].max - loc_v2->coord[1], 2);
    dist = sqrt(dist1) + sqrt(dist2);
    maxDist = dist > maxDist ? dist : maxDist;

    dist1 = pow(MBR[0].max - loc_v1->coord[0], 2) + pow(MBR[1].min - loc_v1->coord[1], 2);
    dist2 = pow(MBR[0].max - loc_v2->coord[0], 2) + pow(MBR[1].min - loc_v2->coord[1], 2);
    dist = sqrt(dist1) + sqrt(dist2);
    maxDist = dist > maxDist ? dist : maxDist;

    dist1 = pow(MBR[0].min - loc_v1->coord[0], 2) + pow(MBR[1].min - loc_v1->coord[1], 2);
    dist2 = pow(MBR[0].min - loc_v2->coord[0], 2) + pow(MBR[1].min - loc_v2->coord[1], 2);
    dist = sqrt(dist1) + sqrt(dist2);
    maxDist = dist > maxDist ? dist : maxDist;

    dist1 = pow(MBR[0].min - loc_v1->coord[0], 2) + pow(MBR[1].max - loc_v1->coord[1], 2);
    dist2 = pow(MBR[0].min - loc_v2->coord[0], 2) + pow(MBR[1].max - loc_v2->coord[1], 2);
    dist = sqrt(dist1) + sqrt(dist2);
    maxDist = dist > maxDist ? dist : maxDist;

    return maxDist;
}

/*
 *	Check whether an object @obj_v is "relevant" to the query @q.
 *	That is, whether @obj_v contains a keyword in the query @q.
 *  modified from is_relevant_obj( obj_v, q)
 */
int is_relevant_obj(obj_t* obj_v, psi_t* psi_v)
{
	int cnt = 0;
	k_node_t* k_node_v;

	k_node_v = psi_v->k_head->next;
	while (k_node_v != NULL) {
		if (has_key_obj(obj_v, k_node_v->key))
			cnt++;

		k_node_v = k_node_v->next;
	}

	return cnt;
}

/*
 *	Release the binary search tree T.
 * modified from bst_release
 */
void release_IF(bst_t* T)
{
	if (T != NULL) {
		if (T->root != NULL)
			release_IF_sub(T->root);

		free(T);

		/*s*/
		stat_v.memory_v -= sizeof(bst_t);
		/*s*/
	}
}

void release_IF_sub(bst_node_t* x)
{
	if (x->left != NULL)
		bst_release_sub(x->left);
	if (x->right != NULL)
		bst_release_sub(x->right);

	release_obj_set(x->p_list_obj);
	free(x);

	/*s*/
	stat_v.memory_v -= sizeof(bst_node_t);
	/*s*/
}

/*
 *	Check the distance constraint.
 */
bool check_dist_constraint(obj_set_t* obj_set_v, obj_t* obj_v, B_KEY_TYPE d)
{
	obj_node_t* obj_node_iter;

	obj_node_iter = obj_set_v->head->next;
	while (obj_node_iter != NULL) {
		if (calc_dist_obj(obj_node_iter->obj_v, obj_v) > d + EPSILON){
			return false;
		}
		obj_node_iter = obj_node_iter->next;
	}

	return true;
}


/*
 *	Check the distance constraint.
 */
bool check_pairwise_dist_constraint(obj_set_t* obj_set_v, B_KEY_TYPE d)
{
	obj_node_t *obj_node_iter, *obj_node_iter2;

	obj_node_iter = obj_set_v->head->next;
	while (obj_node_iter != NULL) {
		obj_node_iter2 = obj_node_iter->next;
		while (obj_node_iter2 != NULL) {
			if (calc_dist_obj(obj_node_iter->obj_v, obj_node_iter2->obj_v) > d + EPSILON)
				return false;
			obj_node_iter2 = obj_node_iter2->next;
		}
		obj_node_iter = obj_node_iter->next;
	}

	return true;
}


/* check whether obj_v is dominated by obj_v2 */
bool isDominated(obj_t* obj_v, obj_t* obj_v2, psi_t* psi_v)
{

	k_node_t* k_node_v;
	B_KEY_TYPE key;

	k_node_v = psi_v->k_head->next;

	//for each keyword in  psi_v
	while (k_node_v != NULL) {
		key = k_node_v->key;
		//if there exist a keyword that is covered by obj_v but not obj_v2
		//obj_v is not dominated by obj_v2
		if (has_key_obj(obj_v, key) && !has_key_obj(obj_v2, key))
			return false;
		k_node_v = k_node_v->next;
	}
	return true;
}


/******
 *	Count the number of keywords that occur in @k_head2 and @k_head1.
 *  modified from "key_exclusion"
 ******/
/* note that there MUST be NO dulipcate in the list
 */
int number_intersection(k_node_t* k_head1, k_node_t* k_head2)
{
	int count = 0;
	k_node_t *k_node_v1, *k_node_v2;
	KEY_TYPE key;

	k_node_v1 = k_head1->next;
	while (k_node_v1 != NULL) {
		key = k_node_v1->key;
		k_node_v2 = k_head2->next;
		while (k_node_v2 != NULL) {
			if (k_node_v2->key == key) {
				count++;
				break;
			}
			k_node_v2 = k_node_v2->next;
		}
		k_node_v1 = k_node_v1->next;
	}

	return count;
}


/*
 *	Exclude the keywords that occur in @obj_v from @psi_v1.
 *
 *	Return the resulting keywords in @psi_v1.
 */
psi_t* psi_exclusion(psi_t* psi_v1, obj_t* obj_v)
{

	psi_t* psi_excluded = alloc_psi();
	k_node_t* k_node_excluded = psi_excluded->k_head; //pointer always point to the end of link list

	int tag;
	k_node_t *k_node_prev, *k_node_v1, *k_node_v2;

	k_node_prev = psi_v1->k_head;
	k_node_v1 = psi_v1->k_head->next;
	while (k_node_v1 != NULL) {
		tag = 0;
		k_node_v2 = obj_v->k_head->next;
		while (k_node_v2 != NULL) {
			if (k_node_v2->key == k_node_v1->key) {
				tag = 1;
				break;
			}

			k_node_v2 = k_node_v2->next;
		}
		if (tag == 1) {
			//The current keyword should be deleted from psi_v1

			//--
			add_keyword_entry(k_node_excluded, k_node_v1->key); //pointed moved to new end inside the function
			psi_excluded->key_n++;
			//--

			k_node_prev->next = k_node_v1->next;
			k_node_v1->next = NULL;
			free(k_node_v1);
			psi_v1->key_n--;

			k_node_v1 = k_node_prev->next;
		} else {
			k_node_v1 = k_node_v1->next;
			k_node_prev = k_node_prev->next;
		}
	}

	return psi_excluded;
}

/*
 */
void psi_insert(psi_t* psi_v, k_node_t* k_head)
{
	k_node_t* k_node_v;

	k_node_v = k_head->next;
	while (k_node_v != NULL) {
		add_keyword_entry(psi_v->k_head, k_node_v->key);
		psi_v->key_n++;

		k_node_v = k_node_v->next;
	}

	return;
}
