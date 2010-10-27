/*
 GTfold: compute minimum free energy of RNA secondary structure
 Copyright (C) 2008  David A. Bader
 http://www.cc.gatech.edu/~bader

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _MAIN_H
#define _MAIN_H

#include <string>
#include <vector>

#include "constants.h"

using namespace std;

//void traverse_wholeData(int[], int, int);

enum GTFOLD_FLAGS { SUCCESS = 0, FAILURE, ERR_OPEN_FILE, NO_CONS_FOUND};

GTFOLD_FLAGS initialize_constraints(int*** fbp, int*** pbp, int& numpConstraints, int& numfConstraints, const char* constr_file);

GTFOLD_FLAGS handle_IUPAC_code(const std::string& s, const int bases);
void limit_contact_distance(int lCD, int length);
void force_noncanonical_basepair(const char* nc_basepairs, int length);

bool is_valid_base(char c)
{	
	return ( (c-'A' == 0) || (c-'a' == 0) || 
			 (c-'C' == 0) || (c-'c' == 0) ||
			 (c-'G' == 0) || (c-'g' == 0) ||
			 (c-'U' == 0) || (c-'u' == 0));
}

void trim_spaces(std::string& str)
{
	// Trim Both leading and trailing spaces
	size_t startpos = str.find_first_not_of(" \t"); // Find the first character position after excluding leading blank spaces
	size_t endpos = str.find_last_not_of(" \t"); // Find the first character position from reverse af

	// if all spaces or empty return an empty string
	if(( string::npos == startpos ) || ( string::npos == endpos))
	{
		str = "";
	}
	else
		str = str.substr( startpos, endpos-startpos+1 );

}

void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters = " ")
{
	// Skip delimiters at beginning.
	std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

	while (std::string::npos != pos || std::string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
}
	
	


#endif
