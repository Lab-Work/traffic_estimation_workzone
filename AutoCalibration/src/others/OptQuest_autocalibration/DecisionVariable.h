#pragma once

class DecisionVariable
{
public:
	DecisionVariable(void);
	~DecisionVariable(void);

	CString m_name;
	long m_ID;
	// Value returned from GetSolutionVariableValue()
	double m_currentValue;  

	// TODO: Add information that will allow you to set values in the simulation
	// model for the entity that corresponds to this decsion variable.  For example, 
	// an Excel cell reference.
};
