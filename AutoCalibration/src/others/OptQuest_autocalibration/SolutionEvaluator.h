#pragma once
#include "ProblemSetup.h"

class SolutionEvaluator
{
public:
	SolutionEvaluator(void);
	~SolutionEvaluator(void);

	bool EvaluateOptQuestSolution(ProblemSetup* setup);

};
