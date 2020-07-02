#include "query1_NonParallel.h"
#include "Incremental_Detector.h"

#include <iostream>
#include <vector>
#include <deque>
#include <string>
#include <array>
#include <map>
#include <thread>
#include <mutex>
#include <chrono>
#include <future>
#include <list>
#include <condition_variable>
#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
#include <curl/curl.h>
#include <unistd.h>

using namespace rapidjson;
using namespace std;

size_t writeFunction_Parallel(void *ptr, size_t size, size_t nmemb, std::string* data) {
    data->append((char*) ptr, size * nmemb);
    return size * nmemb;
}

struct Tuple{
    int sequenceNumber;
    double voltage;
    double current;
};

struct Result{
    int s;
    int d;
    int event_s;
};

int NETWORK_FREQUENCY_Query1 = 50;
int SAMPLERATE_Query1 = 50000;
int period_Query1 = SAMPLERATE_Query1 / NETWORK_FREQUENCY_Query1;
int MAXIMUM_WINDOW_SIZE_Query1 = 100;
extern vector<int> F1Test;

Query1_NonParallel::Query1_NonParallel(string address) {

    cout << "Query 1 starting" << endl;

    Init_dict init_dict;
    Incremental_Detector incremental_Detector(init_dict);

    CURL *curl_get;
    curl_get = curl_easy_init();
    curl_easy_setopt(curl_get, CURLOPT_URL, address.c_str());

    CURL *curl_post;
    curl_post = curl_easy_init();
    curl_easy_setopt(curl_post, CURLOPT_URL, address.c_str());

    while(true){
        std::string response_string;
        curl_easy_setopt(curl_get, CURLOPT_WRITEFUNCTION, writeFunction_Parallel);
        curl_easy_setopt(curl_get, CURLOPT_WRITEDATA, &response_string);
        curl_easy_perform(curl_get);

        if(response_string.empty()){
            usleep(1000000);
            continue;
        }

        if(response_string.size() == 16) {
            cout<< response_string << endl;
            continue;
        }


        Document document;
        document.Parse(response_string.c_str());
        Value::ConstMemberIterator itr = document.FindMember("records");
        if (itr == document.MemberEnd()) {
            cout << response_string << endl;
            usleep(1000000);
            break;
        }

        double voltage[1000];
        double current[1000];
        int batchNumber = itr->value[0]["i"].GetInt() / 1000;
        for (int i = 0; i < itr->value.Size(); i++) {
            voltage[i] = itr->value[i]["voltage"].GetDouble();
            current[i] = itr->value[i]["current"].GetDouble();
        }
        Feature X = incremental_Detector.compute_input_signal(voltage, current, period_Query1, true);
        X.sequenceNumber = batchNumber;

        Event event_interval_indices = incremental_Detector.predict(batchNumber, X);

        Result result_temp;
        result_temp.s = X.sequenceNumber;
        if (event_interval_indices.event_end != 0) {  // if an event is returned
            cout << "batchNumber " <<  X.sequenceNumber << "  Event Detected at   " << event_interval_indices.event_start << "," <<
                 event_interval_indices.event_end << endl;
            int mean_event_index = (event_interval_indices.event_start + event_interval_indices.event_end) / 2 + 1;

            result_temp.d = 1;
            result_temp.event_s = mean_event_index;
            F1Test.push_back(1);

        } else {
            result_temp.d = 0;
            result_temp.event_s = -1;
            incremental_Detector.refresh_window(MAXIMUM_WINDOW_SIZE_Query1);
            F1Test.push_back(0);
        }
        std::string jsonTemp =   "{\"s\":";
        jsonTemp.append(to_string (result_temp.s));
        if(result_temp.d == 0){
            jsonTemp.append(",\"d\":0,\"event_s\":");
        }else{
            jsonTemp.append(",\"d\":1,\"event_s\":");
        }
        if(result_temp.event_s == -1){
            jsonTemp.append(to_string(-1));
        }else {
            jsonTemp.append(to_string(result_temp.event_s));
        }
        jsonTemp.append("}");

        std::string response_string_post;
        struct curl_slist *headers = NULL;
        headers = curl_slist_append(headers, "Content-Type: application/json");
        curl_easy_setopt(curl_post, CURLOPT_HTTPHEADER, headers);
        curl_easy_setopt(curl_post, CURLOPT_POSTFIELDS, jsonTemp.c_str());
        curl_easy_setopt(curl_post, CURLOPT_WRITEFUNCTION, writeFunction_Parallel);
        curl_easy_setopt(curl_post, CURLOPT_WRITEDATA, &response_string_post);
        curl_easy_perform(curl_post);

    }
    cout << "> Query 1 done!" << endl;
}