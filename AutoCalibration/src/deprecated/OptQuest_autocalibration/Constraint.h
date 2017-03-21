#pragma once

class Constraint
{
public:
	Constraint(void);
	~Constraint(void);

	// This class is used to display results at the end of the optimization.
	CString m_Name;
	long m_ID;
};
