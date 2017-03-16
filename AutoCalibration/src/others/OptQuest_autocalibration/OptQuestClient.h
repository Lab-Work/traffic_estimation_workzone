#pragma once
#include "DecisionVariable.h"
#include "SolutionEvaluator.h"
#include "ProblemSetup.h"


class COptQuestClient : public IOptQuestEvents
{
private:
    IOptQuestPtr m_server;
	bool IsOptQuestError();
	long COptQuestClient::GetError();
	void ThrowException();

public: 
	COptQuestClient(void);
	~COptQuestClient(void);

	bool Init();
    // Optimization setup.  Methods to add decision variables, constraints, requirements, and an objective
	long AddContinuousVariable(CString name, double lowerBound, double upperBound);
	long AddDiscreteVariable(CString name, double lowerBound, double upperBound, double stepSize);
	long AddBinaryVariable(CString name);
	long AddIntegerVariable(CString name, double lowerBound, double upperBound);
	long AddDesignVariable(CString name, double lowerBound, double upperBound, double stepSize);
	long AddPermutationVariable(CString name);
	long AddUserControlledVariable(CString name);

	long AddLowerRequirement(CString name, double lowerBound);
	long AddUpperRequirement(CString name, double upperBound);
    long AddDualRequirement(CString name, double lowerBound, double upperBound);
 
 	long AddStringConstraint(CString name, CString expression);
	long AddGEConstraint(CString name, double rhs);
    long AddEQConstraint(CString name, double rhs);
    long AddLEConstraint(CString name, double rhs);
    void AddVariableToConstraint(long constraintID, long variableID, double coeff);

	bool StringConstraintHasUserVariables(long stringID);
	bool StringConstraintIsLinear(long stringID);

	long AddUserControlledObjective(CString name, bool setMaximize);
	long AddStringObjective(CString name, CString expression, bool setMaximize);
 
    // Solution methods.  Methods to set values in a solution and methods to 
	// retrieve information about a solution.
    void SetSolutionVariableValue(long solutionID, long variableID, double value);
	double GetSolutionVariableValue(long solutionID, long variableID);
    void SetSolutionObjectiveValue(long solutionID, long objID, double objValue);
	double GetSolutionObjectiveValue(long solutionID, long objID);
    void SetSolutionRequirementValue(long solutionID, long reqID, double value);
    double GetSolutionRequirementValue(long solutionID, long reqID);
    double GetSolutionConstraintLHSValue(long solutionID, long stringConstraintID);
    double GetSolutionConstraintRHSValue(long solutionID, long stringConstraintID);
    long GetSolutionIteration(long solutionID);
    bool IsSolutionFeasible(long solutionID);
    bool IsSolutionStringConstraintFeasible(long solutionID, long stringConstraintID);
    bool IsSolutionRequirementFeasible(long solutionID, long reqID);


	long GetBestSolution();
	long GetNthBestSolution(long nth);
	long GetTerminationReason();

	// Optimization control
 	void SetLicenseID(long license);
	void SetMaximumIterations(long numIterations);
	void Optimize();
	void StopOptimization();
	void SetMaximumTime(long numSeconds);
	void ContinueOptimize();
	void SetAdditionalIterations(long additonalIterations);
	void SetAdditionalTime(long additionalSeconds);
	void OptimizeAdditional();
    void SetUserControlledStop(bool onOff);
    long GetElapsedTime();
    void SetDatabaseSize(long numSolutions);
    void SetRandomSeed(long seed);
    void LogSetup(CString logFilePath);
    void LogSolutions(CString logFilePath);
    void UseEnglishOnly(bool onOff);
    long GetNumberofCompletedIterations();
    long GetNumberOfFeasibleSolutions();
    long GetNumberOfInfeasibleSolutions();
    void SetCannotGenerateAsTermReason(bool onOff);
    long GetIterationSolution(long iteration);
    void SetUseInfeasibilityIndex(bool onOff);
	void CheckOptimization();

    // EfficientFrontier methods
    long GetNumberEfficientFrontier();
    int GetNthEfficientFrontier(long nth);
    bool IsNthEfficientFrontierFeasible(long nth);

	// Suggested Solutions
	long CreateSolution();
	void AddSuggestedSolution(long solutionID);



    HRESULT __stdcall QueryInterface(const IID &,void **);
    ULONG __stdcall AddRef(void) { return 1; }
    ULONG __stdcall Release(void) { return 1; }
    HRESULT __stdcall EvaluateEvent(long solutionID);
	HRESULT __stdcall MonitorStatusEvent(long solutionID);
	HRESULT __stdcall EfficientFrontierEvent(long solutionID);

	// ProblemSetup is a class that holds the ID's of the objects that were added to the optimization problem
	ProblemSetup* setup;
	SolutionEvaluator evaluator;

	// When there is a com failure, retrieve the IErrorInfo.
	// If it is an OptQuest error, errNumber is an OptQuest 
	// Exception.  See documentation on COptQuestException for the meaning
	// of each error.
	// If it isn't and OptQuest error, errDescription will contrain the exception text
	BSTR errDescription;
	long errNumber;


};
