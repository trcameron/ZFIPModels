#ifndef TESTING_H
#define TESTING_H
#include <stdio.h>
/*
	char_len: character length for pipe calls to nauty
*/
const unsigned int CHAR_LEN = 100;
/*
	geng_call: returns FILE pointer to nauty geng call with given order
*/
FILE* geng_call(const int order,const char *opt);
/*
	genrang_call: returns pipe from nauty genrang call with given order
*/
FILE* genrang_call(const int n,const int num,const int p);
/*
	gentreeg_call: returns pipe from nauty gentreeg call with given order
*/
FILE* gentreeg_call(const int order,const char *opt);
/*
	sizeof_char: number of characters in character array before new line
*/
int sizeof_char(const char* line);
#endif