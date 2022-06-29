#ifndef C_LOADTXT_H
#define C_LOADTXT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
using namespace std;

//------------------------------------------------

template< typename T > vector< vector<T> > loadtxt( const string &filename )
{
   vector< vector<T> > data;
   ifstream in( filename );
   for ( string line; getline( in, line ); )
   {
      stringstream ss( line );
      vector<T> row;
      for ( T d; ss >> d; ) row.push_back( d );
      data.push_back( row );
   }
   return data;
}

//------------------------------------------------

template< typename T > void print( const vector< vector<T> > &data )
{
   for ( auto &row : data )
   {
      for ( auto &item : row ) cout << setw( 10 ) << item << ' ';
      cout << '\n';
   }
}

#endif