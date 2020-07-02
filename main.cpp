#include <iostream>
#include <string>
#include "query1_NonParallel.h"
#include "query2_NonParallel.h"

using namespace std;

int main(int argc, char* argv[]) {

    string QUERY1_ENDPOINT = "/data/1/";
    string QUERY2_ENDPOINT = "/data/2/";

    Query1_NonParallel query1_NonParallel("http://grader/data/1/");
    Query2_NonParallel query2_NonParallel("http://grader/data/2/");
    cout << "Solution Done" << endl;

    return 0;
}