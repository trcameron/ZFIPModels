#include "testing.hpp"
#include <fstream>
#include <iostream>
using namespace std;
/*
	geng_call: returns FILE pointer to nauty geng call with given order
*/
FILE* geng_call(const int order,const char *opt){
	char nauty_call[CHAR_LEN];
	snprintf(nauty_call,CHAR_LEN,"/opt/homebrew/Cellar/nauty/2.9.0/bin/geng %d %s",order,opt);
	return popen(nauty_call, "r");
}
/*
	genrang_call: returns pipe from nauty genrang call with given order
*/
FILE* genrang_call(const int n,const int num,const int p){
	char nauty_call[CHAR_LEN];
	snprintf(nauty_call,CHAR_LEN,"/opt/homebrew/Cellar/nauty/2.9.0/bin/genrang %d %d -g -P%d/10",n,num,p);
	return popen(nauty_call, "r");
}
/*
	gentreeg_call: returns pipe from nauty gentreeg call with given order
*/
FILE* gentreeg_call(const int order,const char *opt){
	char nauty_call[CHAR_LEN];
	snprintf(nauty_call,CHAR_LEN,"/opt/homebrew/Cellar/nauty/2.9.0/bin/gentreeg %d %s",order,opt);
	return popen(nauty_call, "r");
}
/*
	sizeof_char: number of characters in character array before new line
*/
int sizeof_char(const char* line){
	int ind = 0;
	while(line[ind] != '\n'){
		ind++;
	}
	return ind;
}
/*
	main function
*/
/*
int main(){
	// file stream variables
	int status;
	FILE* pipe;
	char line[PATH_MAX];
	string str;
	// call geng for given order
	pipe = geng_call(5,"-g");
	while(fgets(line, PATH_MAX, pipe) != NULL){
		// graph6 string
		for(unsigned int i=0; i<sizeof_char(line); ++i){
			str += line[i];
		}
		cout << str << endl;
		// clear graph6 string
		str.clear();
	}
	// close pipe
	status = pclose(pipe);
	return 0;
}
*/