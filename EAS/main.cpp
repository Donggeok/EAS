#include "eas.h"
#include <Windows.h>

int main() {
	EAS test(100, 50);
	DWORD before_time = GetTickCount();
	test.evolution(EAS_CHOOSE_STRATEGY_CUSTOM, EAS_CROSSOVER_STRATEGY_DOUBLEPOINT, 1.0, EAS_MUTATION_STRATEGY_COMMON, 0.01, 500);
	DWORD after_time = GetTickCount();
	printf("spend time %dms.\n", after_time - before_time);
}