#include <iostream>
#include "eas.h"
#include <ctime>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include "plot.h"

#define M_PI 3.14159265358979323846

using namespace std;


double fitness(__int64 val) {
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

// �Ӵ�С����
bool cmp_fitness(__int64 a, __int64 b) {
	return fitness(a) > fitness(b);
}


class EAS {

private:
	int population_num;
	int coding_len;
	vector<__int64> individuals;

private:
	// �ú������÷������������bit_numλ��С��64λ���������
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
			}
			break;

		case EAS_INIT_STRATEGY_UNIFORM:

			break;

		default:
			break;
		}
	}
	
	void choose_population(int choose_strategy, const vector<__int64> &father, vector<__int64> &child) {

		// һ��Ҫ��ѡ��1/2����Ⱥ
		vector<__int64> tmp_father = father;
		double choose_rate = EAS_CHOOSE_RATE;
		vector<double> fractions(tmp_father.size(), 0.0);

		// ѡ�����
		double sum;
		switch (choose_strategy) {
		case EAS_CHOOSE_STRATEGY_CUSTOM:
			// �Ƚ��о�Ӣ����
			sort(tmp_father.begin(), tmp_father.end(), cmp_fitness);
			child.insert(child.end(), tmp_father.begin(), tmp_father.begin() + EAS_CHOOSE_RESERVE_RATE*population_num);
			// �������̶�(Ӧ�ð�������ľ�Ӣ)
			sum = 0.0;
			
			for (int i = 0; i < tmp_father.size(); ++i) {
				sum += fitness(tmp_father[i]);
				fractions[i] = sum;
			}
			for (int i = 0; i < tmp_father.size(); ++i) {
				fractions[i] = fractions[i] / sum;
			}
			// ɸѡ
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
			// ɸѡ
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
				// 2i��2i+1�ĸ������н���
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

	void print_in_format(__int64 val, double *values) {
		__int64 x1, x2;
		int x1len = 18, x2len = 15;
		double x1a = -3.0, x1b = 12.1, x2a = 4.1, x2b = 5.8;
		x1 = val >> x2len;
		x2 = (x1 << x2len) ^ val;
		double x1real = x1a + (x1b - x1a)*x1 / ((1LL << x1len) - 1LL);
		double x2real = x2a + (x2b - x2a)*x2 / ((1LL << x2len) - 1LL);
		double result = 21.5 + x1real * sin(4 * M_PI * x1real) + x2real * sin(20 * M_PI * x2real);
		// ������ֵ����ȥ
		*values = result;
		printf("x1: %lf, x2: %lf, f(x1, x2): %lf\n", x1real, x2real, result);
	}

public:
	
	// ������ʱ����Ѿ�����Ⱥ��ʼ�����
	EAS(int n, int len, int strategy = EAS_INIT_STRATEGY_RANDOM) {
		population_num = n;
		coding_len = len;
		individuals.resize(population_num);

		init_population(strategy);
	}

	// Ϊ����Ӧ���������Ӧ�Ⱥ��������Խ�����Ϊһ������ʹ��
	// ����ʱ����Ҫ��ѡ����ԣ�������ԣ�������ʣ�������ԣ�������ʣ�
	void evolution(int choose_strategy, int crossover_strategy, double crossover_rate,
					int mutation_strategy, double mutation_rate, int times) {

		vector<__int64> father_before, child, all, results;
		double *best_values = new double[times];
		for (int i = 0; i < times; ++i) {
			// ��ʼ��
			father_before.clear();
			child.clear();
			all.clear();
			results.clear();

			// ѡ��
			choose_population(EAS_CHOOSE_STRATEGY_ROULETTE, individuals, father_before);
			// ����
			crossover_population(crossover_strategy, crossover_rate, father_before, child);
			// ����
			mutation_population(mutation_strategy, mutation_rate, child);
			// ѡ���Ӵ�
			all = individuals;
			all.insert(all.end(), child.begin(), child.end());
			choose_population(choose_strategy, all, results);
			individuals.clear();
			individuals = results;

			//sort(individuals.begin(), individuals.end(), cmp_fitness);
			print_in_format(individuals[0], &best_values[i]);
			
		}
		// ��ͼ
		cvNamedWindow("����ֵ������仯ͼ", 1);
		CPlot plot;
		plot.x_max = (double)times; //�����趨��������������Сֵ
		plot.x_min = 0.0;
		plot.y_max = 39.0;
		plot.y_min = 36.0;
		plot.axis_color = Scalar(255, 0, 0);
		plot.text_color = Scalar(255, 255, 255);
		plot.plot(best_values, times, CV_RGB(0, 0, 0)); //����ֻ����Yֵ XĬ�ϴ�0��ʼ 
		plot.title("f(x1, x2)"); //�����趨���� ֻ����Ӣ�� ���Ļ����� �н�����������Ǻ��鷳
		plot.xlabel("generation", Scalar(255, 255, 0));
		plot.ylabel("f(x1, x2)", Scalar(255, 255, 0));
		cvShowImage("����ֵ������仯ͼ", plot.Figure);
		cvWaitKey(0);
		delete[] best_values;
	}
};

int main() {
	EAS test(300, 33);
	test.evolution(EAS_CHOOSE_STRATEGY_CUSTOM, EAS_CROSSOVER_STRATEGY_DOUBLEPOINT, 0.9, EAS_MUTATION_STRATEGY_COMMON, 0.01, 50);
}