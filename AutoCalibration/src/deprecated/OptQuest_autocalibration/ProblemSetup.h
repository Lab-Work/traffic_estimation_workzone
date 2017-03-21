#pragma once

#include "DecisionVariable.h"
#include "Constraint.h"

class ProblemSetup
{
public:
	ProblemSetup(void);
	~ProblemSetup(void);

	// TODO: The following lists are used to track the decision variables, constraints, requirements and objectives
	// that have been added to the optimization problem.   These should be defined in a way that is meaningful to your
	// particular implementation.
	CArray<DecisionVariable, DecisionVariable&> m_decVars;
	CArray<Constraint, Constraint&> m_constraints;

	long m_objID;
	double m_objValue;

};
