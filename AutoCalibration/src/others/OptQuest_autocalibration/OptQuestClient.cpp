#include "StdAfx.h"
#include "OptQuestClient.h"
#include <iostream>


COptQuestClient::COptQuestClient(void)
{
}

COptQuestClient::~COptQuestClient(void)
{
}

/*
 * Tell the server you are listening for events.  
 */
bool COptQuestClient::Init()
{
    HRESULT hr = S_OK;

    hr = m_server.CreateInstance(__uuidof(OptQuestServer));
    ATLASSERT(SUCCEEDED(hr));
	bool success = true;
   
	hr = m_server->AddListener(this);
	if (SUCCEEDED(hr))
		return true;
	else
		return false;
}

HRESULT COptQuestClient::QueryInterface(const IID & iid,void ** pp)
{
    if (iid == __uuidof(IOptQuestEvents) ||
        iid == __uuidof(IUnknown))
    {
        *pp = this;
        AddRef();
        return S_OK;
    }
    return E_NOINTERFACE;
}

/*
 *  EvaluateEvent() is called when OptQuest has a solution ready for evaluation.
 *	
 */
HRESULT COptQuestClient::EvaluateEvent(long solutionID)
{
	double varValue;
	bool error = false;
	HRESULT hr;
	// Get the decision variable values
	for (int i = 0; i < setup->m_decVars.GetSize(); i++)
	{
		hr = m_server->GetSolutionVariableValue(solutionID, setup->m_decVars.GetAt(i).m_ID, &varValue);
		if (SUCCEEDED(hr))
			setup->m_decVars.GetAt(i).m_currentValue = varValue;
		else
		{
			error = true;
			break;
		}
	}
	bool success = evaluator.EvaluateOptQuestSolution(setup);


	if (success)
	{
		// Set the objective value in the solution
		hr = m_server->SetSolutionObjectiveValue(solutionID, setup->m_objID, setup->m_objValue);
		if (SUCCEEDED(hr))
		{
			hr = m_server->EvaluateComplete(solutionID);
			if (SUCCEEDED(hr))
				return S_OK;
			else
				error = true;
		}		
		else
			error = true;
	}

	if (!success)
	{
		// TODO: The string description will have the error information.  
		// It is probably a COptQuestException indicating a bad ID was used in one of the calls.
		IErrorInfo* pErrInfo = NULL;
		HRESULT hr2 = GetErrorInfo(0,&pErrInfo);
		BSTR description;
		hr2 = pErrInfo->GetDescription(&description);
		return E_FAIL;
	}
    return S_OK;
}

/*
 * MonitorStatusEvent() is called when OptQuest has completed evaluating the solution.  This
 * means the solution has been checked for feasibility and has been checked to see if it is
 * a new, best solution.  If you are providing status information to the user at each iteration
 * you should do it here.  
 */
HRESULT COptQuestClient::MonitorStatusEvent(long solutionID)
{
	// TODO: Provide status information to the user such as the number of completed iterations
	// and the best solution found thus far.
	printf("\nOptQuest finished iteration %d, with best solution: \n", GetNumberofCompletedIterations());

	// get the best solution so far
	// long best_solutionID = m_server->getBestSolution()

	return S_OK;
}

/*
 * EfficientFrontierEvent() is called when OptQuest has completed evaluating solutions for a point on 
 * the Efficient Frontier.  The solutionID is the ID of the best solution found thus far for that point.
 * The method GetNumberEfficientFrontier() returns the number of points that have been evaluated.
 * The method GetNthEfficientFrontier() allows you to retrieve the best solution for each point on the
 * Efficient Frontier. As the optimization continues it it possible that better solutions are found for
 * previously evaluted points on the Efficient Frontier.
 * The method IsNthEfficientFrontierFeasible() allows you to find out if a feasible solution was found
 * for each point on the Efficient Frontier
 */
HRESULT COptQuestClient::EfficientFrontierEvent(long solutionID)
{
	// TODO: Update EfficientFrontier information
	return S_OK;
}

/*
 * Checks the optimization setup for errors which are returned as exceptions.
 * User can call StringConstraintIsLinear() method and StringConstraintHasUserConrolledVariables()
 * after calling this method.
 */
void COptQuestClient::CheckOptimization()
{
	HRESULT hr = m_server->CheckOptimization();
	if (SUCCEEDED(hr))
		return;

	// FAILED 
	ThrowException();

}
long COptQuestClient::AddContinuousVariable(CString name, double lowerBound, double upperBound)
{
	long ID = -1;
	BSTR bstrName = name.AllocSysString();

	HRESULT hr = m_server->AddContinuousVariable(bstrName, lowerBound, upperBound, &ID);
	::SysFreeString(bstrName);

	if (SUCCEEDED(hr))
		return ID;

	// FAILED 
	ThrowException();
	return ID;
}

long COptQuestClient::AddDiscreteVariable(CString name, double lowerBound, double upperBound, double stepSize)
{
	long ID = -1;
	BSTR bstrName = name.AllocSysString();

	HRESULT hr = m_server->AddDiscreteVariable(bstrName, lowerBound, upperBound, stepSize, &ID);
	::SysFreeString(bstrName);

	if (SUCCEEDED(hr))
		return ID;
	// FAILED
	ThrowException();
	return ID;
}

long COptQuestClient::AddBinaryVariable(CString name)
{
	long ID = -1;
	BSTR bstrName = name.AllocSysString();

	HRESULT hr = m_server->AddBinaryVariable(bstrName, &ID);
	::SysFreeString(bstrName);

	if (SUCCEEDED(hr))
		return ID;

	// FAILED 
	ThrowException();
	return ID;
}

long COptQuestClient::AddIntegerVariable(CString name, double lowerBound, double upperBound)
{
	long ID = -1;
	BSTR bstrName = name.AllocSysString();

	HRESULT hr = m_server->AddIntegerVariable(bstrName, lowerBound, upperBound, &ID);
	::SysFreeString(bstrName);

	if (SUCCEEDED(hr))
		return ID;
	// FAILED
	ThrowException();
	return ID;
}

long COptQuestClient::AddDesignVariable(CString name, double lowerBound, double upperBound, double stepSize)
{
	long ID = -1;
	BSTR bstrName = name.AllocSysString();

	HRESULT hr = m_server->AddDesignVariable(bstrName, lowerBound, upperBound, stepSize, &ID);
	::SysFreeString(bstrName);

	if (SUCCEEDED(hr))
		return ID;
	// FAILED
	ThrowException();
	return ID;
}

long COptQuestClient::AddPermutationVariable(CString name)
{
	long ID = -1;
	BSTR bstrName = name.AllocSysString();

	HRESULT hr = m_server->AddPermutationVariable(bstrName,&ID);
	::SysFreeString(bstrName);

	if (SUCCEEDED(hr))
		return ID;
	// FAILED
	ThrowException();
	return ID;
}

long COptQuestClient::AddUserControlledVariable(CString name)
{
	long ID = -1;
	BSTR bstrName = name.AllocSysString();

	HRESULT hr = m_server->AddUserControlledVariable(bstrName, &ID);
	::SysFreeString(bstrName);

	if (SUCCEEDED(hr))
		return ID;
	// FAILED
	ThrowException();
	return ID;
}

long COptQuestClient::AddStringConstraint(CString name, CString expression)
{
	long ID = -1;
	BSTR bstrName = name.AllocSysString();
	BSTR bstrExp = expression.AllocSysString();

	HRESULT hr = m_server->AddStringConstraint(bstrName, bstrExp, &ID);
	::SysFreeString(bstrName);
	::SysFreeString(bstrExp);

	if (SUCCEEDED(hr))
		return ID;
	// FAILED
	ThrowException();
	return ID;
}

bool COptQuestClient::StringConstraintHasUserVariables(long stringConstraintID)
{
	VARIANT_BOOL yesNo;
	HRESULT hr = m_server->StringConstraintHasUserVariables(stringConstraintID, &yesNo);
	if (SUCCEEDED(hr))
		return yesNo == VARIANT_TRUE;
	// FAILED
	ThrowException();
	return false;
}

bool COptQuestClient::StringConstraintIsLinear(long stringConstraintID)
{
	VARIANT_BOOL yesNo;
	HRESULT hr = m_server->StringConstraintIsLinear(stringConstraintID, &yesNo);
	if (SUCCEEDED(hr))
		return yesNo == VARIANT_TRUE;
	// FAILED
	ThrowException();
	return false;
}

long COptQuestClient::AddLowerRequirement(CString name, double lowerBound)
{
	long ID = -1;
	BSTR bstrName = name.AllocSysString();

	HRESULT hr = m_server->AddLowerRequirement(bstrName, lowerBound, &ID);
	::SysFreeString(bstrName);

	if (SUCCEEDED(hr))
		return ID;
	// FAILED
	ThrowException();
	return ID;
}
long COptQuestClient::AddUpperRequirement(CString name, double upperBound)
{
	long ID = -1;
	BSTR bstrName = name.AllocSysString();

	HRESULT hr = m_server->AddUpperRequirement(bstrName, upperBound, &ID);
	::SysFreeString(bstrName);

	if (SUCCEEDED(hr))
		return ID;
	// FAILED
	ThrowException();
	return ID;
}

long COptQuestClient::AddDualRequirement(CString name, double lowerBound, double upperBound)
{
	long ID = -1;
	BSTR bstrName = name.AllocSysString();

	HRESULT hr = m_server->AddDualRequirement(bstrName, lowerBound, upperBound, &ID);
	::SysFreeString(bstrName);

	if (SUCCEEDED(hr))
		return ID;
	// FAILED
	ThrowException();
	return ID;
}

long COptQuestClient::AddGEConstraint(CString name, double rhs)
{
	long ID = -1;
	BSTR bstrName = name.AllocSysString();

	HRESULT hr = m_server->AddGEConstraint(bstrName, rhs, &ID);
	::SysFreeString(bstrName);

	if (SUCCEEDED(hr))
		return ID;
	// FAILED
	ThrowException();
	return ID;
}

long COptQuestClient::AddEQConstraint(CString name, double rhs)
{
	long ID = -1;
	BSTR bstrName = name.AllocSysString();

	HRESULT hr = m_server->AddEQConstraint(bstrName, rhs, &ID);
	::SysFreeString(bstrName);

	if (SUCCEEDED(hr))
		return ID;
	// FAILED
	ThrowException();
	return ID;
}

long COptQuestClient::AddLEConstraint(CString name, double rhs)
{
	long ID = -1;
	BSTR bstrName = name.AllocSysString();

	HRESULT hr = m_server->AddLEConstraint(bstrName, rhs, &ID);
	::SysFreeString(bstrName);

	if (SUCCEEDED(hr))
		return ID;
	// FAILED
	ThrowException();
	return ID;
}

void COptQuestClient::AddVariableToConstraint(long constraintID, long variableID, double coeff)
{
	HRESULT hr = m_server->AddVariableToConstraint(constraintID, variableID, coeff);

	if (SUCCEEDED(hr))
		return;
	// FAILED
	ThrowException();
}


void COptQuestClient::SetLicenseID(long license)
{
	HRESULT hr = m_server->SetLicenseID(license);
	if (SUCCEEDED(hr))
		return;
	// FAILED 
	ThrowException();
}

long COptQuestClient::AddUserControlledObjective(CString name, bool setMaximize)
{
	long ID = -1;
	BSTR bstrName = name.AllocSysString();
	VARIANT_BOOL vSetMaximize;
	if (setMaximize)
		vSetMaximize = VARIANT_TRUE;
	else
		vSetMaximize = VARIANT_FALSE;

	HRESULT hr = m_server->AddUserControlledObjective(bstrName, vSetMaximize, &ID);
	::SysFreeString(bstrName);

	if (SUCCEEDED(hr))
		return ID;
	// FAILED
	ThrowException();
	return ID;
}

long COptQuestClient::AddStringObjective(CString name, CString expression, bool setMaximize)
{
	long ID = -1;
	BSTR bstrName = name.AllocSysString();
	BSTR bstrExp = expression.AllocSysString();
	VARIANT_BOOL vSetMaximize;
	if (setMaximize)
		vSetMaximize = VARIANT_TRUE;
	else
		vSetMaximize = VARIANT_FALSE;

	HRESULT hr = m_server->AddStringObjective(bstrName, bstrExp, vSetMaximize, &ID);
	::SysFreeString(bstrName);
	::SysFreeString(bstrExp);

	if (SUCCEEDED(hr))
		return ID;
	// FAILED
	ThrowException();
	return ID;
}

void COptQuestClient::SetMaximumIterations(long numIterations)
{
	HRESULT hr = m_server->SetMaximumIterations(numIterations);
	if (SUCCEEDED(hr))
		return;
	// FAILED 
	ThrowException();
}

void COptQuestClient::Optimize()
{
	HRESULT hr = m_server->Optimize();
	if (SUCCEEDED(hr))
		return;
	// FAILED 
	ThrowException();
}

void COptQuestClient::SetSolutionVariableValue(long solutionID, long variableID, double value)
{
	HRESULT hr = m_server->SetSolutionVariableValue(solutionID, variableID, value);
	if (SUCCEEDED(hr))
		return;
	// FAILED
	ThrowException();
}

double COptQuestClient::GetSolutionVariableValue(long solutionID, long variableID)
{
	double varValue;
	HRESULT hr = m_server->GetSolutionVariableValue(solutionID, variableID, &varValue);
	if (SUCCEEDED(hr))
		return varValue;
	// FAILED
	ThrowException();
	return 0.0;
}

/*
 *  Set the value of the objective in the solution.  
 */
void COptQuestClient::SetSolutionObjectiveValue(long solutionID, long objID, double objValue)
{
	HRESULT hr = m_server->SetSolutionObjectiveValue(solutionID, objID, objValue);
	if (SUCCEEDED(hr))
		return;
	// FAILED 
	ThrowException();
}

/*
 * Get the value of the objective in the solution
 */
double COptQuestClient::GetSolutionObjectiveValue(long solutionID, long objID)
{
	double objValue;
	HRESULT hr = m_server->GetSolutionObjectiveValue(solutionID, objID, &objValue);
	if (SUCCEEDED(hr))
		return objValue;
	// FAILED
	ThrowException();
	return 0.0;
}

/*
 * Set the value of the requirement in the solution
 */
void COptQuestClient::SetSolutionRequirementValue(long solutionID, long reqID, double value)
{
	HRESULT hr = m_server->SetSolutionRequirementValue(solutionID, reqID, value);
	if (SUCCEEDED(hr))
		return;
	// FAILED 
	ThrowException();
}
/*
 * Get the value of the requirement in the solution.
 */
double COptQuestClient::GetSolutionRequirementValue(long solutionID, long reqID)
{
	double reqValue;
	HRESULT hr = m_server->GetSolutionRequirementValue(solutionID, reqID, &reqValue);
	if (SUCCEEDED(hr))
		return reqValue;
	// FAILED
	ThrowException();
	return 0.0;
}

/*
 * Calculates the value of the left-hand-side of the string constraint using the values in the solution.
 */
double COptQuestClient::GetSolutionConstraintLHSValue(long solutionID, long stringConstraintID)
{
	double lhsValue;
	HRESULT hr = m_server->GetSolutionConstraintLHSValue(solutionID, stringConstraintID, &lhsValue);
	if (SUCCEEDED(hr))
		return lhsValue;
	// FAILED
	ThrowException();
	return 0.0;
}
/*
 * Calculates the value of the right-hand-side of the string constraint using the values in the solution.
 */
double COptQuestClient::GetSolutionConstraintRHSValue(long solutionID, long stringConstraintID)
{
	double rhsValue;
	HRESULT hr = m_server->GetSolutionConstraintRHSValue(solutionID, stringConstraintID, &rhsValue);
	if (SUCCEEDED(hr))
		return rhsValue;
	// FAILED
	ThrowException();
	return 0.0;
}
/*
 * Returns the iteration of the solution.  Note: the iteration number is being used as the solutionID
 * so as long as the solutionID is greater than 0, it can be used as the iteration number.
 */
long COptQuestClient::GetSolutionIteration(long solutionID)
{
	return solutionID;
	return -1;
}
/*
 * Returns true if the solution satisfies all constraints and requirements
 */
bool COptQuestClient::IsSolutionFeasible(long solutionID)
{
	VARIANT_BOOL feasible;
	HRESULT hr = m_server->IsSolutionFeasible(solutionID, &feasible);
	if (SUCCEEDED(hr))
		return feasible == VARIANT_TRUE;
	// FAILED
	ThrowException();
	return false;
}
bool COptQuestClient::IsSolutionStringConstraintFeasible(long solutionID, long stringConstraintID)
{
	VARIANT_BOOL feasible;
	HRESULT hr = m_server->IsSolutionStringConstraintFeasible(solutionID, stringConstraintID, &feasible);
	if (SUCCEEDED(hr))
		return feasible == VARIANT_TRUE;
	// FAILED
	ThrowException();
	return false;
}
bool COptQuestClient::IsSolutionRequirementFeasible(long solutionID, long reqID)
{
	VARIANT_BOOL feasible;
	HRESULT hr = m_server->IsSolutionRequirementFeasible(solutionID, reqID, &feasible);
	if (SUCCEEDED(hr))
		return feasible == VARIANT_TRUE;
	// FAILED
	ThrowException();
	return false;
}
/*
 * Get the solutionID of the best solution found thus far.  
 */
long COptQuestClient::GetBestSolution()
{
	long solutionID;
	HRESULT hr = m_server->GetBestSolution(&solutionID);
	if (SUCCEEDED(hr))
		return solutionID;
	// FAILED
	ThrowException();
	return -1;
}

/* Returns the reason the optimization stopped.  
 * 0 = The optimization has not started
 * 1 = The optimization is still running	 
 * 3 = The optimization was solved using an Linear/Integer/Mixed Integer Program.
 * 4 = The optimization stopped due to the Auto Stop feature
 * 5 = All solutions have been generated.
 * 6 = The optimization stopped when the maximum number of iterations was reached.
 * 7 = The optimization stopped when the maximum time was reached.
 * 8 = The optimization was stopped by the user.
 * 10 = The optimization stopped due to an exception.
 * 11 = The Meta-Heuristic factories have completed.
 * 12 = There are no solutions that satisfy the constraints.
 * 13 = New (different) solutions cannot be generated.
 */

long COptQuestClient::GetTerminationReason()
{
	long reason;
	HRESULT hr = m_server->GetTerminationReason(&reason);
	if (SUCCEEDED(hr))
		return reason;
	// FAILED!!
	ThrowException();
	return -1;
}
/*
 *
 * Stops the optimization.  If the current solution has not completed evaluation, it is
 * put back onto the queue of unevaluated solutions.  It will be the first solution evaluated
 * if the optimization is continued.
*/
void COptQuestClient::StopOptimization()
{
	HRESULT hr = m_server->StopOptimization();
	if (SUCCEEDED(hr))
		return ;
	// FAILED!!
	ThrowException();
}
/*
 * Sets the maximum number of seconds the optimization should run.  The user can stop
 * the optimization before the time has elapsed by calling StopOptimization()
 */
void COptQuestClient::SetMaximumTime(long numSeconds)
{
	HRESULT hr = m_server->SetMaximumTime(numSeconds);
	if (SUCCEEDED(hr))
		return ;
	// FAILED!!
	ThrowException();
}
/*
 * Continues an optimization that was stopped before the maximum time or maximum iterations 
 * had been reached.
 */
void COptQuestClient::ContinueOptimize()
{
	HRESULT hr = m_server->ContinueOptimize();
	if (SUCCEEDED(hr))
		return ;
	// FAILED!!
	ThrowException();
}

/*
 * An optimization that has completed can be continued by calling this method and then 
 * OptimizeAdditonal() to run the additional iterations
 */
void COptQuestClient::SetAdditionalIterations(long additonalIterations)
{
	HRESULT hr = m_server->SetAdditionalIterations(additonalIterations);
	if (SUCCEEDED(hr))
		return ;
	// FAILED!!
	ThrowException();
}
/*
 * An optimization that has completed can be continued by calling this method and then 
 * OptimizeAdditonal() to run for the additional seconds
 */
void COptQuestClient::SetAdditionalTime(long additionalSeconds)
{
	HRESULT hr = m_server->SetAdditionalTime(additionalSeconds);
	if (SUCCEEDED(hr))
		return ;
	// FAILED!!
	ThrowException();
}
/*
 * This method allows the user to run additional iterations or seconds after an optimization has completed.
 * To run additional iterations, call the method SetAdditionalIterations() specifying the number of additional iterations to be run.  
 * To run additional seconds, call the method SetAdditionalTime() specifying the number of additional seconds to be run.
 */
void COptQuestClient::OptimizeAdditional()
{
	HRESULT hr = m_server->OptimizeAdditional();
	if (SUCCEEDED(hr))
		return ;
	// FAILED!!
	ThrowException();
}

/*
 * Returns the solution ID of the nth best solution where solutions are ordered from
 * best to worst.  If the input parameter is 1, the best solution is returned,  If the input 
 * parameter is 2, the second best solution etc. If the nth best 
 * solution doesn't exist, an exception is thrown.
 */
long COptQuestClient::GetNthBestSolution(long nth)
{
	long solutionID;
	HRESULT hr = m_server->GetNthBestSolution(nth, &solutionID);
	if (SUCCEEDED(hr))
		return solutionID;
	// FAILED!!
	ThrowException();
	return -1;
}

/*
 * If the input parameter is true, the optimization is going to be stopped by calling StopOptimization()
 * rather than by setting a number of iterations (SetMaximumIterations()) or a maximum number of seconds 
  * to run (SetMaximumTime()).
 */
void COptQuestClient::SetUserControlledStop(bool onOff)
{
	VARIANT_BOOL vOnOff;
	if (onOff)
		vOnOff = VARIANT_TRUE;
	else
		vOnOff = VARIANT_FALSE;
	HRESULT hr = m_server->SetUserControlledStop(vOnOff);
	if (SUCCEEDED(hr))
		return;
	// FAILED!!
	ThrowException();
}

/*
 * Returns the number of seconds the optimization has been running.
 */
long COptQuestClient::GetElapsedTime()
{
	long elapsedSeconds;
	HRESULT hr = m_server->GetElapsedTime(&elapsedSeconds);
	if (SUCCEEDED(hr))
		return elapsedSeconds;
	// FAILED!!
	ThrowException();
	return -1;
}
/*
 * Sets the number of solutions OptQuest should keep in its internal database of solutions.
 */
void COptQuestClient::SetDatabaseSize(long numSolutions)
{
	HRESULT hr = m_server->SetDatabaseSize(numSolutions);
	if (SUCCEEDED(hr))
		return;
	// FAILED!!
	ThrowException();
}

/*
 * OptQuest uses a random seed in solution generation.  This allows you
 * to set the random seed. 
 */
void COptQuestClient::SetRandomSeed(long seed)
{
	HRESULT hr = m_server->SetDatabaseSize(seed);
	if (SUCCEEDED(hr))
		return;
	// FAILED!!
	ThrowException();
}
/*
 * When this method is called, OptQuest logs the optimization setup as an xml file.
 * Often used to debug setup problems.
 */
void COptQuestClient::LogSetup(CString logFilePath)
{
	BSTR bstrPath = logFilePath.AllocSysString();

	HRESULT hr = m_server->LogSetup(bstrPath);
	::SysFreeString(bstrPath);

	if (SUCCEEDED(hr))
		return;

	// FAILED 
	ThrowException();
}

/*
 * Logs all solutions in a csv formatted file.  This is used for debugging solution generation.
 */
void COptQuestClient::LogSolutions(CString logFilePath)
{
	BSTR bstrPath = logFilePath.AllocSysString();

	HRESULT hr = m_server->LogSolutions(bstrPath);
	::SysFreeString(bstrPath);

	if (SUCCEEDED(hr))
		return;

	// FAILED 
	ThrowException();
}

/*
 * When set to true OptQuest does not use localization to parse numbers used in string
 * constraints or string objectives.  All numbers are treated as being enterd in English.
 * Excel has an English only setting (somewhere).  If you set this in Excel, you need 
 * to use this method to notify OptQuest. 
 */
void COptQuestClient::UseEnglishOnly(bool onOff)
{
	
	VARIANT_BOOL vOnOff;
	if (onOff)
		vOnOff = VARIANT_TRUE;
	else
		vOnOff = VARIANT_FALSE;
	HRESULT hr = m_server->UseEnglishOnly(vOnOff);
	if (SUCCEEDED(hr))
		return;
	// FAILED!!
	ThrowException();
}

/*
 * Returns the number of iterations OptQuest has completed.
 */
long COptQuestClient::GetNumberofCompletedIterations()
{
	long iterations;
	HRESULT hr = m_server->GetNumberofCompletedIterations(&iterations);
	if (SUCCEEDED(hr))
		return iterations;
	// FAILED!!
	ThrowException();
	return -1;
}
/*
 * Returns the number of solutions that satisfied all constraints and requirements.
 */
long COptQuestClient::GetNumberOfFeasibleSolutions()
{
	long numFeasible;
	HRESULT hr = m_server->GetNumberOfFeasibleSolutions(&numFeasible);
	if (SUCCEEDED(hr))
		return numFeasible;
	// FAILED!!
	ThrowException();
	return -1;
}
/*
 * Returns the number of solutions that violated one or more constraints or requirements.
 */
long COptQuestClient::GetNumberOfInfeasibleSolutions()
{
	long numInfeasible;
	HRESULT hr = m_server->GetNumberOfInfeasibleSolutions(&numInfeasible);
	if (SUCCEEDED(hr))
		return numInfeasible;
	// FAILED!!
	ThrowException();
	return -1;
}

/*
 * The default behavior when OptQuest cannot generate anymore solutions is to throw exception -14000.
 * Setting the input parameter to true overrides this behavior and OptQuest returns termination reason
 * 13 rather than throwing an exception.
 */
void COptQuestClient::SetCannotGenerateAsTermReason(bool onOff)
{
	VARIANT_BOOL vOnOff;
	if (onOff)
		vOnOff = VARIANT_TRUE;
	else
		vOnOff = VARIANT_FALSE;
	HRESULT hr = m_server->SetCannotGenerateAsTermReason(vOnOff);
	if (SUCCEEDED(hr))
		return;
	// FAILED!!
	ThrowException();
}

/*
 * Returns the solution ID for the solution evaluated at the specified iteration.
 */
long COptQuestClient::GetIterationSolution(long iteration)
{
	long solution;
	HRESULT hr = m_server->GetIterationSolution(iteration, &solution);
	if (SUCCEEDED(hr))
		return solution;
	// FAILED!!
	ThrowException();
	return -1;
}

/*
 * When this feature is turned on, solutions are ranked according to their proximity to feasibility which affects the order 
 * solutions are returned when GetNthBestSolution() is called.
 */
void COptQuestClient::SetUseInfeasibilityIndex(bool onOff)
{
	VARIANT_BOOL vOnOff;
	if (onOff)
		vOnOff = VARIANT_TRUE;
	else
		vOnOff = VARIANT_FALSE;
	HRESULT hr = m_server->SetUseInfeasibilityIndex(vOnOff);
	if (SUCCEEDED(hr))
		return;
	// FAILED!!
	ThrowException();
}

/*
 * Returns the number of Efficient Frontier points that have been evaluated.
 */
long COptQuestClient::GetNumberEfficientFrontier()
{
	long count;
	HRESULT hr = m_server->GetNumberEfficientFrontier(&count);
	if (SUCCEEDED(hr))
		return count;
	// FAILED!!
	ThrowException();
	return -1;
}

/*
 * Returns the solution ID of the best solution for the nth Efficient Frontier
 * point where nth begins at 1.
 */
int COptQuestClient::GetNthEfficientFrontier(long nth)
{
	long solution;
	HRESULT hr = m_server->GetNthEfficientFrontier(nth, &solution);
	if (SUCCEEDED(hr))
		return solution;
	// FAILED!!
	ThrowException();
	return -1;
}

/*
 * Creates a new solution that can be used to create a suggested solution. The variable values are 
 * initialized to the last suggested value or the midpoint of the lower and upper bound if no solution has been suggested. 
 * Use SetSolutionVariableValue() to set the values of each variable.  Call AddSuggestedSolution() when the variables have
 * been set.
 */
long COptQuestClient::CreateSolution()
{
	long solutionID;
	HRESULT hr = m_server->CreateSolution(&solutionID);
	if (SUCCEEDED(hr))
		return solutionID;
	// FAILED!!
	ThrowException();
	return -1;
}
/*
 * The input solution is added to the set of suggested solutions. If the optimization has not started, the 
 * suggested solution will be one of the first solutions evaluated. If the optimization is running, the 
 * suggested solution will be one of the next solutions to be evaluated. 
 * After this method has been called, the solutionID is no longer valid and should not be used to access
 * the solution.
 */
void COptQuestClient::AddSuggestedSolution(long solutionID)
{
	HRESULT hr = m_server->AddSuggestedSolution(solutionID);
	if (SUCCEEDED(hr))
		return;
	// FAILED!!
	ThrowException();
}

/*
 * Returns true if the best solution for the nth point of the Efficient Frontier
 * satisfied all constraints and requirements. The parameter nth begins at 1.
 */
bool COptQuestClient::IsNthEfficientFrontierFeasible(long nth)
{	VARIANT_BOOL feasible;
	HRESULT hr = m_server->IsNthEfficientFrontierFeasible(nth, &feasible);
	if (SUCCEEDED(hr))
		return feasible == VARIANT_TRUE;
	// FAILED!!
	ThrowException();
	return false;
}

bool COptQuestClient::IsOptQuestError()
{
	// OptQuest errors are numeric values between -1200 and -16104.
	// If converting to an integer fails, or the value is not within
	// this range it is not an OptQuest exception
	IErrorInfo* pErrInfo = NULL;
	HRESULT hr = GetErrorInfo(0,&pErrInfo);
	hr = pErrInfo->GetDescription(&errDescription);
	HRESULT hr2 = VarI4FromStr(errDescription, 0,0, &errNumber);
	if (SUCCEEDED(hr2))
	{
		if (errNumber >= -16104 && errNumber <= -1200)
			return true;
	}
	return false;
}
/*
 * GetError() is called when HRESULT indicates an error. 
 * Use the IErrorInfo interface to get the description
 * which is the error number from the COptQuestException.
 * TODO: Make this smarter.  If the description is not a
 * COptQuestException then we want the description.
 */
//long COptQuestClient::GetError()
//{
//	IErrorInfo* pErrInfo = NULL;
//	HRESULT hr = GetErrorInfo(0,&pErrInfo);
//	hr = pErrInfo->GetDescription(&errDescription);
//	HRESULT hr2 = VarI4FromStr(description, 0,0, &errNumber);
//	if (SUCCEEDED(hr2))
//		return err;
//	else
//		return hr2;
//}

void COptQuestClient::ThrowException()
{
	if (IsOptQuestError())
		throw(errNumber);
	else
		throw(errDescription);
}


