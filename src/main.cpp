/*
 *	Author: Harry Kai-Ho Chan
 *	Email: kai-ho@ruc.dk
 */


#include <iostream>
#include <fstream>
#include "cdcoskq.h"

using namespace std;

IRTree_t IRTree_v;
coskq_stat_t stat_v;
bst_t* IF_global;

int cost_tag;
int w_tag;
B_KEY_TYPE MAX_DIA;
W_TYPE MAX_WEIGHT;

///
data_t* data_v;
bool firstTime = true;
unordered_map<KEY_TYPE, KEY_TYPE>* keyfreq_hashmap;
//memory.
float first_memory_v;
float first_memory_max;

float first_tree_memory_v;
float first_tree_memory_max;

//bool* nullquery;
void cdcoskq(bool prune_tag);

void build_IF(data_t* data_v);

int main(int argc, char* argv[])
{
	cdcoskq(true);
    return 0;
}

void cdcoskq(bool prune_tag)
{
    int i, size_sum;
    B_KEY_TYPE cost, cost_sum, weight_sum, B_sum, ratioa_sum, ratiob_sum, ratiob2_sum;
	W_TYPE weight;
	coskq_config_t* cfg;
    //coskq_stat_t* sta_v;
    query_t** q_set;
    range* MBR;
    obj_set_t *S, *S_o;
    FILE* r_fp;

    //--

    memset(&stat_v, 0, sizeof(coskq_stat_t));

    size_sum = 0;
    //Read the cofig.
    ///read the config.txt, return coskq_config_t pointer
    printf("Reading configuration ...\n");
    cfg = read_config_coskq();

    cfg->prune_opt = 0;

    cost_tag = cfg->cost_measure;
    w_tag = cfg->w_opt;

    printf(" --- ");

    if (cost_tag == 1)
        printf("MaxSum\t");
    else if (cost_tag == 2)
        printf("Dia\t");
    else {
        printf("cost function undefined..");
        return;
    }

    if (w_tag == 1)
        printf("MaxWeight\t");
    else if (w_tag == 2)
        printf("SumWeight\t");
    else {
        printf("weight function undefined..");
        return;
    }

	if (cfg->alg_opt == 1){
        printf("CD-Exact:\t");
		cfg->m_opt = 1;
	}else if (cfg->alg_opt == 2){
        printf("CD-Appro:\t");
		cfg->m_opt = 2;
	}else if (cfg->alg_opt == 3){
		printf("Appro Adapt:\t");
		if (cfg->m_opt == 1)
			printf("Cao-Appro\t");
		else if (cfg->m_opt == 2)
			printf("Long-Appro\t");
		else {
			printf("method option undefined..");
			return;
		}
	}else if (cfg->alg_opt == 4)
        printf("Baseline Exact:\t");
//    else if (cfg->alg_opt == 5)
//        printf("Baseline Appro:\t");
//    else if (cfg->alg_opt == 7)
//		printf("Appro Adapt:\t");
//	else if (cfg->alg_opt == 8)
//		printf("Long-Exact:\t");
    else {
        printf("algorithm option undefined..");
        return;
    }

    printf(" --- \n");

    //Read the data.
    ///read the "loc" and "doc" files, return data_t pointer
    if (firstTime) {
        printf("Reading data ...\n");
        keyfreq_hashmap = new unordered_map<KEY_TYPE, KEY_TYPE>();
		MAX_WEIGHT = 0;
        data_v = read_data_coskq(cfg, keyfreq_hashmap, MAX_WEIGHT);
//		printf("MAX_WEIGHT:%d\n", MAX_WEIGHT);
//		printf("#obj:%d\n", data_v->obj_n);
//        printf("#key (in hashmap):%lu\n", keyfreq_hashmap->size());
    }
#ifndef WIN32
    float sys_t, usr_t, usr_t_sum = 0;
    struct rusage IR_tree_sta, IR_tree_end;

    GetCurTime(&IR_tree_sta);
#endif
  
    if (firstTime) {

        printf("Building IR-tree ...\n");
        build_IRTree(data_v);
        //IR-tree is pointed by global var IRTree_v now
        //print_and_check_tree( 1, cfg->tree_file);
        //check_IF( );
        first_memory_v = stat_v.memory_v;
        first_memory_max = stat_v.memory_max;
        first_tree_memory_v = stat_v.tree_memory_v;
        first_tree_memory_max = stat_v.tree_memory_max;
    } else {
        stat_v.memory_v += first_memory_v;
        stat_v.memory_max += first_memory_max;
        stat_v.tree_memory_v += first_tree_memory_v;
        stat_v.tree_memory_max += first_tree_memory_max;
    }

#ifndef WIN32
    GetCurTime(&IR_tree_end);
    GetTime(&IR_tree_sta, &IR_tree_end, &stat_v.irtree_build_time, &sys_t);
#endif

    if (firstTime) {
        firstTime = false;
        //Build Inverted file
        build_IF(data_v);
    }
    //Get the whole range.
    MBR = get_MBR_node(IRTree_v.root, IRTree_v.dim);

    //---
    double x = MBR[0].max - MBR[0].min;
    double y = MBR[1].max - MBR[1].min;
    double maxD = sqrt(x * x + y * y);
//    printf("maxD:%0.3lf\n", maxD);
	MAX_DIA = maxD;
    //---

    //Generate the set of querys.
    ///within the MBR
    printf("Generating queries ...\n");

	if(cfg->low == 0)
		q_set = gen_query_set_Cao(cfg->q_set_size, cfg->q_key_n, cfg->obj_n, data_v);
	else
		q_set = gen_query_set2(cfg->q_set_size, cfg->q_key_n, MBR, data_v, cfg->low, cfg->high);
	
    if (q_set == NULL) {
        printf("Query generation failed!\n");
        exit(0);
    }
	
    //Query.
    printf("Performing Queries ...\n");

    
    cost_sum = 0;
    weight_sum = 0;
    B_sum = 0;
    ratioa_sum = 0;
    ratiob_sum = 0;
    ratiob2_sum = 0;

    int numOfSnotNull = 0;
    int numOfExecutedQuery = 0;
    int numOfExactNotNull = 0;
	
    ///for each query
    for (i = 0; i < cfg->q_set_size; i++) {

        printf("Query #%i ...\n", i + 1);

        double maxQ = calc_maxDist(MBR, q_set[i]->loc_v);
		B_KEY_TYPE B ;
		B_KEY_TYPE cost_S_c =-1;
//		if(cfg->B >99)
		{
			obj_set_t* S_c = CostEnum_Appro(q_set[i]);
			cost_S_c = comp_cost(cost_tag, S_c, q_set[i]);
			B = cost_S_c* ((cfg->B) / 100.0);
			release_obj_set(S_c);
		}
		
        numOfExecutedQuery++;

#ifndef WIN32
        struct rusage query_sta, query_end;

        GetCurTime(&query_sta);
#endif
		
        if (cfg->alg_opt == 1)
            S = Exact(q_set[i], cfg->m_opt, B);
        else if (cfg->alg_opt == 2)
            S = Appro(q_set[i], cfg->m_opt, B);
        else if (cfg->alg_opt == 3)
			S = ApproAdapt(q_set[i],cfg->m_opt,B);
        else if (cfg->alg_opt == 4)
            S = ExactBaseline(q_set[i], cfg->m_opt, B);
//        else if (cfg->alg_opt == 5)
//			S = NULL;
//		 else if (cfg->alg_opt == 6)
//		S= NULL;
//		else if (cfg->alg_opt==7)
//			S = ApproAdapt(q_set[i],cfg->m_opt,B);
//		else if(cfg->alg_opt==8) //Long-Exact
//			S = CostEnum(q_set[i], cfg->m_opt, 0);
//
        //debug purpose.
//		  S_o = Exact(q_set[i], cfg->m_opt, B);
//                if (S != NULL && S_o != NULL && comp_weight(w_tag, S) != comp_weight(w_tag, S_o))
//        	    		{
//                			print_obj_set(S_o, stdout);
//                			printf("weight:%d\tcost:%.3lf\n", comp_weight(w_tag, S_o),comp_cost(cost_tag,S_o, q_set[i]));
//
//                			print_obj_set(S, stdout);
//                        			printf("weight:%d\tcost:%.3lf\n", comp_weight(w_tag, S),comp_cost(cost_tag,S, q_set[i]));
//                                }

        if (S == NULL) {
            //printf("\n S = NULL\n");
            S = alloc_obj_set(); //avoid null pointer later
            //            return;

        } else {
            numOfSnotNull++; //number of divide for average weight, cost
        }

#ifndef WIN32
        GetCurTime(&query_end);

        GetTime(&query_sta, &query_end, &usr_t, &sys_t);
        usr_t_sum += usr_t;
#endif
		cost = comp_cost(cost_tag, S, q_set[i]);
        cost_sum += cost;
		
		weight = comp_weight(w_tag, S);
		
		weight_sum += weight;
        size_sum += S->obj_n;
        B_sum += B;

        //-----------
        //Ratio processing.
        B_KEY_TYPE weight_opt = 0, cost_opt = 0;
        B_KEY_TYPE ratioB2;

		if(prune_tag &&( cfg->alg_opt==2||cfg->alg_opt==3)){//only appro
        obj_set_t* S2 = Exact(q_set[i], 1, B); //a ratio
//			obj_set_t* S2 = S_o;
        if (S2 != NULL) {
            B_KEY_TYPE ratioA = 0, ratioB = 0;
			
			cost_opt = comp_cost(cost_tag, S2, q_set[i]);
			weight_opt = comp_weight(w_tag, S2);
			
			
            if (weight_opt != 0)
                ratioA = (float)weight / (float)weight_opt;
            if (cost_opt != 0)
                ratioB = cost / cost_opt;

            ratioa_sum += (float)ratioA;
            ratiob_sum += (float)ratioB;

            numOfExactNotNull++; //number to divide for average a ratio
        }
        release_obj_set(S2);
		}

        ratioB2 = cost / B;
        ratiob2_sum += (float)ratioB2;
        //--------

        //Print the accumulated query results.
        if (i == 0) {
            if ((r_fp = fopen(COSKQ_RESULT_FILE, "w")) == NULL) {
                fprintf(stderr, "Cannot open the coskq_result file.\n");
                exit(0);
            }

        } else {
            if ((r_fp = fopen(COSKQ_RESULT_FILE, "a")) == NULL) {
                fprintf(stderr, "Cannot open the coskq_result file.\n");
                exit(0);
            }
        }

        fprintf(r_fp, "Query #%i:\n", i + 1);
        fprintf(r_fp, "Keywords: ");
        print_k_list(q_set[i]->psi_v->k_head, r_fp);
        ///
        fprintf(r_fp, "Locations: ");
        print_loc(q_set[i]->loc_v, r_fp);
        ///

        //Print the query result.
        if (S->obj_n > 0)
            print_obj_set(S, r_fp);
        else
            fprintf(r_fp, "\nS = NULL\n\n");

        fprintf(r_fp, "Dist: %0.4lf\nB:%0.4lf\n", cost, B);
        //        printf("Cost: %0.4lf\n", cost);
        printf("Dist: %0.4lf\nB:%0.4lf\n", cost, B);
        fprintf(r_fp, "Weight: %d\n", weight);
        printf("Weight: %d\n", weight);
#ifndef WIN32
        fprintf(r_fp, "Time: %f\n\n", usr_t);
        printf("Time: %f\n\n", usr_t);
#endif

        fclose(r_fp);

        //Print the statistics.
        ///it is updated after each quey performed
        ///k>1 here, no need to k+1 here.
        if (numOfSnotNull != 0) {
            stat_v.aver_cost = cost_sum / (numOfSnotNull);
            stat_v.aver_weight = weight_sum / (numOfSnotNull);
            stat_v.ratiob2_aver = ratiob2_sum / (numOfSnotNull);

        } else {
            stat_v.aver_cost = 0;
            stat_v.aver_weight = 0;
            stat_v.ratiob2_aver = 0;
        }

        if (numOfExactNotNull != 0) {
            stat_v.aver_B = B_sum / (numOfExactNotNull);
            stat_v.ratioa_aver = ratioa_sum / (numOfExactNotNull);
            stat_v.ratiob_aver = ratiob_sum / (numOfExactNotNull);
        } else {
            stat_v.aver_B = 0;
            stat_v.ratioa_aver = 0;
            stat_v.ratiob_aver = 0;
        }

#ifndef WIN32
        if (numOfExecutedQuery != 0)
            stat_v.q_time = usr_t_sum / (numOfExecutedQuery);
        else
            stat_v.q_time = 0;
#endif
        print_coskq_stat(cfg, numOfExecutedQuery);

        release_obj_set(S);

    } ///end for each query

    //Release the resources.
    for (i = 0; i < cfg->q_set_size; i++)
        release_query(q_set[i]);
    free(q_set);
    free(MBR);
    //printf("Memory balance: %f\n", stat_v.memory_v / cfg->q_set_size);
    free(cfg);
}

/*
 * Construct the inverted file @IF_v based on the data
 */
void build_IF(data_t* data_v)
{
//	bst_t* IF_global; //comment this line if this global index is to be used
	
    bst_node_t* bst_node_v;

    IF_global = bst_ini();

    //Insert all the objects to construct the IF
    for (int i = 0; i < data_v->obj_n; i++) {

        k_node_t* k_node_v = data_v->obj_v[i].k_head->next;
        while (k_node_v != NULL) {
            bst_node_v = bst_search(IF_global, k_node_v->key);

            if (bst_node_v != NULL) {
                add_obj_set_entry(&data_v->obj_v[i], bst_node_v->p_list_obj);
            } else //bst_node_v = NULL
            {
                //Insert a new keyword in the binary tree.
                bst_node_v = (bst_node_t*)malloc(sizeof(bst_node_t));
                memset(bst_node_v, 0, sizeof(bst_node_t));

                /*s*/
                stat_v.memory_v += sizeof(bst_node_t);
                if (stat_v.memory_v > stat_v.memory_max)
                    stat_v.memory_max = stat_v.memory_v;
                /*s*/

                //Update the posting list.
                bst_node_v->key = k_node_v->key;
                bst_node_v->p_list_obj = alloc_obj_set();

                add_obj_set_entry(&data_v->obj_v[i], bst_node_v->p_list_obj);
                bst_insert(IF_global, bst_node_v);
            }
            k_node_v = k_node_v->next;
        }
    }
}
