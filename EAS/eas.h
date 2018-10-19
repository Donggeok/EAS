#pragma once

#include <iostream>
#include "eas.h"
#include <ctime>
#include <vector>
#include <algorithm>
#include <cstdlib>

#ifndef CV_SWAP
#define CV_SWAP(a,b,t) ((t) = (a), (a) = (b), (b) = (t))
#endif

#ifndef MIN
#  define MIN(a,b)  ((a) > (b) ? (b) : (a))
#endif

#define CV_IMPLEMENT_QSORT_EX( func_name, T, LT, user_data_type )                   \
void func_name( T *array, size_t total, user_data_type aux )                        \
{                                                                                   \
    int isort_thresh = 7;                                                           \
    T t;                                                                            \
    int sp = 0;                                                                     \
                                                                                    \
    struct                                                                          \
    {                                                                               \
        T *lb;                                                                      \
        T *ub;                                                                      \
    }                                                                               \
    stack[48];                                                                      \
                                                                                    \
    (void)aux;                                                                      \
                                                                                    \
    if( total <= 1 )                                                                \
        return;                                                                     \
                                                                                    \
    stack[0].lb = array;                                                            \
    stack[0].ub = array + (total - 1);                                              \
                                                                                    \
    while( sp >= 0 )                                                                \
    {                                                                               \
        T* left = stack[sp].lb;                                                     \
        T* right = stack[sp--].ub;                                                  \
                                                                                    \
        for(;;)                                                                     \
        {                                                                           \
            int i, n = (int)(right - left) + 1, m;                                  \
            T* ptr;                                                                 \
            T* ptr2;                                                                \
                                                                                    \
            if( n <= isort_thresh )                                                 \
            {                                                                       \
            insert_sort:                                                            \
                for( ptr = left + 1; ptr <= right; ptr++ )                          \
                {                                                                   \
                    for( ptr2 = ptr; ptr2 > left && LT(ptr2[0],ptr2[-1]); ptr2--)   \
                        CV_SWAP( ptr2[0], ptr2[-1], t );                            \
                }                                                                   \
                break;                                                              \
            }                                                                       \
            else                                                                    \
            {                                                                       \
                T* left0;                                                           \
                T* left1;                                                           \
                T* right0;                                                          \
                T* right1;                                                          \
                T* pivot;                                                           \
                T* a;                                                               \
                T* b;                                                               \
                T* c;                                                               \
                int swap_cnt = 0;                                                   \
                                                                                    \
                left0 = left;                                                       \
                right0 = right;                                                     \
                pivot = left + (n/2);                                               \
                                                                                    \
                if( n > 40 )                                                        \
                {                                                                   \
                    int d = n / 8;                                                  \
                    a = left, b = left + d, c = left + 2*d;                         \
                    left = LT(*a, *b) ? (LT(*b, *c) ? b : (LT(*a, *c) ? c : a))     \
                                      : (LT(*c, *b) ? b : (LT(*a, *c) ? a : c));    \
                                                                                    \
                    a = pivot - d, b = pivot, c = pivot + d;                        \
                    pivot = LT(*a, *b) ? (LT(*b, *c) ? b : (LT(*a, *c) ? c : a))    \
                                      : (LT(*c, *b) ? b : (LT(*a, *c) ? a : c));    \
                                                                                    \
                    a = right - 2*d, b = right - d, c = right;                      \
                    right = LT(*a, *b) ? (LT(*b, *c) ? b : (LT(*a, *c) ? c : a))    \
                                      : (LT(*c, *b) ? b : (LT(*a, *c) ? a : c));    \
                }                                                                   \
                                                                                    \
                a = left, b = pivot, c = right;                                     \
                pivot = LT(*a, *b) ? (LT(*b, *c) ? b : (LT(*a, *c) ? c : a))        \
                                   : (LT(*c, *b) ? b : (LT(*a, *c) ? a : c));       \
                if( pivot != left0 )                                                \
                {                                                                   \
                    CV_SWAP( *pivot, *left0, t );                                   \
                    pivot = left0;                                                  \
                }                                                                   \
                left = left1 = left0 + 1;                                           \
                right = right1 = right0;                                            \
                                                                                    \
                for(;;)                                                             \
                {                                                                   \
                    while( left <= right && !LT(*pivot, *left) )                    \
                    {                                                               \
                        if( !LT(*left, *pivot) )                                    \
                        {                                                           \
                            if( left > left1 )                                      \
                                CV_SWAP( *left1, *left, t );                        \
                            swap_cnt = 1;                                           \
                            left1++;                                                \
                        }                                                           \
                        left++;                                                     \
                    }                                                               \
                                                                                    \
                    while( left <= right && !LT(*right, *pivot) )                   \
                    {                                                               \
                        if( !LT(*pivot, *right) )                                   \
                        {                                                           \
                            if( right < right1 )                                    \
                                CV_SWAP( *right1, *right, t );                      \
                            swap_cnt = 1;                                           \
                            right1--;                                               \
                        }                                                           \
                        right--;                                                    \
                    }                                                               \
                                                                                    \
                    if( left > right )                                              \
                        break;                                                      \
                    CV_SWAP( *left, *right, t );                                    \
                    swap_cnt = 1;                                                   \
                    left++;                                                         \
                    right--;                                                        \
                }                                                                   \
                                                                                    \
                if( swap_cnt == 0 )                                                 \
                {                                                                   \
                    left = left0, right = right0;                                   \
                    goto insert_sort;                                               \
                }                                                                   \
                                                                                    \
                n = MIN( (int)(left1 - left0), (int)(left - left1) );               \
                for( i = 0; i < n; i++ )                                            \
                    CV_SWAP( left0[i], left[i-n], t );                              \
                                                                                    \
                n = MIN( (int)(right0 - right1), (int)(right1 - right) );           \
                for( i = 0; i < n; i++ )                                            \
                    CV_SWAP( left[i], right0[i-n+1], t );                           \
                n = (int)(left - left1);                                            \
                m = (int)(right1 - right);                                          \
                if( n > 1 )                                                         \
                {                                                                   \
                    if( m > 1 )                                                     \
                    {                                                               \
                        if( n > m )                                                 \
                        {                                                           \
                            stack[++sp].lb = left0;                                 \
                            stack[sp].ub = left0 + n - 1;                           \
                            left = right0 - m + 1, right = right0;                  \
                        }                                                           \
                        else                                                        \
                        {                                                           \
                            stack[++sp].lb = right0 - m + 1;                        \
                            stack[sp].ub = right0;                                  \
                            left = left0, right = left0 + n - 1;                    \
                        }                                                           \
                    }                                                               \
                    else                                                            \
                        left = left0, right = left0 + n - 1;                        \
                }                                                                   \
                else if( m > 1 )                                                    \
                    left = right0 - m + 1, right = right0;                          \
                else                                                                \
                    break;                                                          \
            }                                                                       \
        }                                                                           \
    }                                                                               \
}

#define hough_cmp_gt(l1,l2) (aux[l1] > aux[l2])

static CV_IMPLEMENT_QSORT_EX(sort_only_index, int, hough_cmp_gt, const float*)

// 初始化策略
#define EAS_INIT_STRATEGY_RANDOM    0
#define EAS_INIT_STRATEGY_UNIFORM   1

// 交叉策略
#define EAS_CROSSOVER_STRATEGY_SINGLEPOINT   0
#define EAS_CROSSOVER_STRATEGY_DOUBLEPOINT   1
#define EAS_CROSSOVER_STRATEGY_MULTIPOINT    2
#define EAS_CROSSOVER_STRATEGY_UNIFORM       3
#define EAS_CROSSOVER_STRATEGY_PARTIALMAPPED 4
#define EAS_CROSSOVER_STRATEGY_ORDER         5
#define EAS_CROSSOVER_STRATEGY_CYCLE         6

// 变异策略
#define EAS_MUTATION_STRATEGY_COMMON     0
#define EAS_MUTATION_STRATEGY_EXCHANGE   1
#define EAS_MUTATION_STRATEGY_SHIFT      2
#define EAS_MUTATION_STRATEGY_INSERT     3
#define EAS_MUTATION_STRATEGY_UNIFORM    4
#define EAS_MUTATION_STRATEGY_NORMAL     5
#define EAS_MUTATION_STRATEGY_NO_UNIFORM 6

// 选择策略
#define EAS_CHOOSE_STRATEGY_ROULETTE     0
#define EAS_CHOOSE_STRATEGY_ORDER        1
#define EAS_CHOOSE_STRATEGY_CROWING      2
#define EAS_CHOOSE_STRATEGY_COMPITIVE    3
#define EAS_CHOOSE_STRATEGY_CUSTOM       4


// 其他参数
// 该参数基本不要动，设置的小了之后父代减少，产生的子代更少，导致种群越来越少（类比与现实生活中的生子率，维持1稳定）
#define EAS_CHOOSE_RATE           1.0
// 该参数代表精英保留预留率
#define EAS_CHOOSE_RESERVE_RATE   0.1


#define M_PI 3.14159265358979323846


#define BAG_ITEM_NUM 50
#define BAG_VOLUME   1000
int C[BAG_ITEM_NUM] = { 220, 208, 198, 192, 180, 180, 165, 162, 160, 158,
			155, 130, 125, 122, 120, 118, 115, 110, 105, 101, 
			100, 100, 98, 96, 95, 90, 88, 82, 80, 77,
			75, 73, 72, 70, 69, 66, 65, 63, 60, 58, 
			56, 50, 30, 20, 15, 10, 8, 5, 3, 1 };
int W[BAG_ITEM_NUM] = { 80, 82, 85, 70, 72, 70, 66, 50, 55, 25,
			50, 55, 40, 48, 50, 32, 22, 60, 30, 32, 
			40, 38, 35, 32, 25, 28, 30, 22, 50, 30, 
			45, 30, 60, 50, 20, 65, 20, 25, 30, 10, 
			20, 25, 15, 10, 10, 10, 4, 4, 2, 1 };

using namespace std;

void greedy_operator(__int64 &val) {
	float b[BAG_ITEM_NUM];
	int b_size = 0;
	vector<int> sort_buf(BAG_ITEM_NUM, 0);
	for (int i = 0; i < BAG_ITEM_NUM; ++i) {
		if (val & (1LL << i)) {
			b[i] = ((float)C[i] / W[i]);
			sort_buf[b_size] = i;
			b_size++;
		}
	}
	sort_only_index(&sort_buf[0], b_size, b);
	int cur_v = 0;
	bool no_enough = true;
	for (int i = 0; i < b_size; ++i) {
		if (no_enough) {
			cur_v += W[sort_buf[i]];
			if (cur_v > BAG_VOLUME) {
				cur_v -= W[sort_buf[i]];
				val &= ~(1LL << sort_buf[i]);
				no_enough = false;
			}
		}
		else {
			val &= ~(1LL << sort_buf[i]);
		}
	}
}

double fitness_raw(__int64 val) {
	__int64 x1, x2;
	int x1len = 18, x2len = 15;
	double x1a = -3.0, x1b = 12.1, x2a = 4.1, x2b = 5.8;
	x1 = val >> x2len;
	x2 = (x1 << x2len) ^ val;
	double x1real = x1a + (x1b - x1a)*x1 / ((1LL << x1len) - 1LL);
	double x2real = x2a + (x2b - x2a)*x2 / ((1LL << x2len) - 1LL);
	double result = 21.5 + x1real * sin(4 * M_PI * x1real) + x2real * sin(20 * M_PI * x2real);
	return result;
}

int fitness(__int64 val) {
	int result = 0;
	for (int i = 0; i < BAG_ITEM_NUM; ++i) {
		if (val & (1LL << i)) {
			result += C[i];
		}
	}
	return result;
}

// 从大到小排序
bool cmp_fitness(__int64 a, __int64 b) {
	return fitness(a) > fitness(b);
}


class EAS {

private:
	int population_num;
	int coding_len;
	vector<__int64> individuals;

private:
	// 该函数的用法就是随机生成bit_num位（小于64位）的随机数
	__int64 Rand64(int bit_num) {
		__int64 nRandData = 0;
		for (int i = 0; i < 8; ++i) {
			int nMask = 256;
			__int64 nData = rand() % nMask;
			nRandData |= (nData << (i * 8));
		}
		__int64 mask = 0LL;
		for (int i = 0; i < bit_num; ++i) {
			mask |= 1LL << i;
		}
		nRandData &= mask;
		return nRandData;
	}

	void init_population(int init_strategy) {
		switch (init_strategy) {
		case EAS_INIT_STRATEGY_RANDOM:
			srand((unsigned int)time(NULL));
			for (int i = 0; i < individuals.size(); ++i) {
				individuals[i] = Rand64(coding_len);
				greedy_operator(individuals[i]);
			}
			break;

		case EAS_INIT_STRATEGY_UNIFORM:

			break;

		default:
			break;
		}
	}

	void choose_population(int choose_strategy, const vector<__int64> &father, vector<__int64> &child) {

		// 一般要求选择1/2的种群
		vector<__int64> tmp_father = father;
		double choose_rate = EAS_CHOOSE_RATE;
		vector<double> fractions(tmp_father.size(), 0.0);
		int t;

		// 选择策略
		double sum;
		switch (choose_strategy) {
		case EAS_CHOOSE_STRATEGY_CUSTOM:
			// 先进行精英保留
			sort(tmp_father.begin(), tmp_father.end(), cmp_fitness);
			child.insert(child.end(), tmp_father.begin(), tmp_father.begin() + EAS_CHOOSE_RESERVE_RATE*population_num);
			// 进行轮盘赌(包括上面的精英)
			sum = 0.0;

			// 记t为轮盘赌起始的位置
			t = 0;

			for (int i = t; i < tmp_father.size(); ++i) {
				sum += fitness(tmp_father[i]);
				fractions[i] = sum;
			}
			for (int i = t; i < tmp_father.size(); ++i) {
				fractions[i] = fractions[i] / sum;
			}
			// 筛选
			srand((unsigned int)time(NULL));
			for (int i = 0; i < (EAS_CHOOSE_RATE - EAS_CHOOSE_RESERVE_RATE)*population_num; ++i) {
				double randDouble = rand() / (RAND_MAX + 1.0);
				int j;
				for (j = 0; j < population_num; ++j) {
					if (randDouble < fractions[j]) {
						break;
					}
				}
				child.push_back(tmp_father[j]);
			}
			break;

		case EAS_CHOOSE_STRATEGY_ROULETTE:
			sum = 0.0;

			for (int i = 0; i < tmp_father.size(); ++i) {
				sum += fitness(tmp_father[i]);
				fractions[i] = sum;
			}
			for (int i = 0; i < tmp_father.size(); ++i) {
				fractions[i] = fractions[i] / sum;
			}
			// 筛选
			srand((unsigned int)time(NULL));
			for (int i = 0; i < EAS_CHOOSE_RATE*population_num; ++i) {
				double randDouble = rand() / (RAND_MAX + 1.0);
				int j;
				for (j = 0; j < population_num; ++j) {
					if (randDouble < fractions[j]) {
						break;
					}
				}
				child.push_back(tmp_father[j]);
			}
			break;
		default:
			break;
		}
	}

	void crossover_population(int crossover_strategy, double crossover_rate,
		const vector<__int64> &father, vector<__int64> &child) {
		int tmp_time = father.size() / 2;
		srand((unsigned int)time(NULL));
		for (int i = 0; i < tmp_time; ++i) {
			double randDouble = rand() / (RAND_MAX + 1.0);
			if (randDouble < crossover_rate) {
				// 2i和2i+1的父辈进行交叉
				int p1, p2;
				__int64 p_mask, n_mask, childa, childb;
				switch (crossover_strategy) {
				case EAS_CROSSOVER_STRATEGY_SINGLEPOINT:

					break;
				case EAS_CROSSOVER_STRATEGY_DOUBLEPOINT:
					p1 = rand() % (coding_len - 1) + 1;
					p2 = rand() % (coding_len - 1) + 1;
					if (p2 == p1) {
						p2 = (p1 + coding_len / 2) % (coding_len - 1) + 1;
					}
					if (p1 > p2) {
						int tmp = p1;
						p1 = p2;
						p2 = tmp;
					}

					p_mask = 0LL;
					for (int j = p1; j < p2; ++j) {
						p_mask |= (1LL << j);
					}
					n_mask = ~p_mask;
					childa = p_mask&father[2 * i] | n_mask&father[2 * i + 1];
					childb = p_mask&father[2 * i + 1] | n_mask&father[2 * i];

					child.push_back(childa);
					child.push_back(childb);

					break;
				default:
					break;
				}
			}
			else {
				continue;
			}
		}
	}

	void mutation_population(int mutation_strategy, double mutation_rate, vector<__int64> &child) {
		switch (mutation_strategy) {
		case EAS_MUTATION_STRATEGY_COMMON:
			srand((unsigned int)time(NULL));
			for (int i = 0; i < child.size(); ++i) {
				for (int j = 0; j < coding_len; ++j) {
					double randDouble = rand() / (RAND_MAX + 1.0);
					if (randDouble < mutation_rate) {
						child[i] ^= (1LL << j);
					}
				}
			}
			break;

		default:
			break;
		}
	}

	void print_in_format_bak(__int64 val, double *values) {
		__int64 x1, x2;
		int x1len = 18, x2len = 15;
		double x1a = -3.0, x1b = 12.1, x2a = 4.1, x2b = 5.8;
		x1 = val >> x2len;
		x2 = (x1 << x2len) ^ val;
		double x1real = x1a + (x1b - x1a)*x1 / ((1LL << x1len) - 1LL);
		double x2real = x2a + (x2b - x2a)*x2 / ((1LL << x2len) - 1LL);
		double result = 21.5 + x1real * sin(4 * M_PI * x1real) + x2real * sin(20 * M_PI * x2real);
		// 将最后的值传出去
		*values = result;
		printf("x1: %.15lf, x2: %.15lf, f(x1, x2): %.15lf\n", x1real, x2real, result);
	}

	void print_in_format(__int64 val) {
		printf("选择的物品编号：\n");
		int count_sum = 0;
		int weight_sum = 0;
		for (int i = 0; i < BAG_ITEM_NUM; ++i) {
			if (val & (1LL << i)) {
				printf("%d ", i + 1);
				count_sum += C[i];
				weight_sum += W[i];
			}
		}
		printf("\n");
		printf("总价值：%d，总重量：%d\n", count_sum, weight_sum);
	}

public:

	// 创建的时候就已经将种群初始化完毕
	EAS(int n, int len, int strategy = EAS_INIT_STRATEGY_RANDOM) {
		population_num = n;
		coding_len = len;
		individuals.resize(population_num);

		init_population(strategy);
	}

	// 为了适应更多的求适应度函数，所以将其作为一个函数使用
	// 进化时还需要有选择策略，交叉策略，交叉概率，变异策略，变异概率，
	void evolution(int choose_strategy, int crossover_strategy, double crossover_rate,
		int mutation_strategy, double mutation_rate, int times) {

		vector<__int64> father_before, child, all, results;
		int i;
		for (i = 0; i < times; ++i) {
			// 初始化
			father_before.clear();
			child.clear();
			all.clear();
			results.clear();

			// 选择
			choose_population(EAS_CHOOSE_STRATEGY_ROULETTE, individuals, father_before);
			// 交叉
			crossover_population(crossover_strategy, crossover_rate, father_before, child);
			// 变异
			mutation_population(mutation_strategy, mutation_rate, child);
			// 选择子代
			all = individuals;
			all.insert(all.end(), child.begin(), child.end());

			for (int i = 0; i < all.size(); ++i) {
				greedy_operator(all[i]);
			}

			choose_population(choose_strategy, all, results);
			individuals.clear();
			individuals = results;

			print_in_format(individuals[0]);
		}
	}
};