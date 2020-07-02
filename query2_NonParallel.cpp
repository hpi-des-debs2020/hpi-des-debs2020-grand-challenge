#include "query2_NonParallel.h"
#include "Incremental_Detector_Disorder.h"

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
#include <condition_variable>
#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
#include <curl/curl.h>
#include <unistd.h>
#include <list>

using namespace rapidjson;
using namespace std;

size_t writeFunction_Parallel2(void *ptr, size_t size, size_t nmemb, std::string* data) {
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

int NETWORK_FREQUENCY_Query2 = 50;
int SAMPLERATE_Query2 = 50000;
int period_Query2 = SAMPLERATE_Query2 / NETWORK_FREQUENCY_Query2;
int MAXIMUM_WINDOW_SIZE_Query2 = 100;

int Radius_Query2 =7;
const int queueSize_Query2 = 30;
int maxThrehold_Query2 = 1000;
int lateTime_Query2 = 20;

Query2_NonParallel::Query2_NonParallel(string address) {
    cout << "Query 2 starting" << endl;

    Init_dict init_dict;
    Incremental_Detector_Disorder incremental_Detector_Disorder(init_dict);
    vector<Tuple> inputBatch_Query2[queueSize_Query2];
    list<int> waitBatch;

    CURL *curl_get;
    curl_get = curl_easy_init();
    curl_easy_setopt(curl_get, CURLOPT_URL, address.c_str());

    CURL *curl_post;
    curl_post = curl_easy_init();
    curl_easy_setopt(curl_post, CURLOPT_URL, address.c_str());

    int batchNumber = 0;

    while(true){
        std::string response_string;
        curl_easy_setopt(curl_get, CURLOPT_WRITEFUNCTION, writeFunction_Parallel2);
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
        for (int i = 0; i < itr->value.Size(); i++) {
            struct Tuple tuple;
            tuple.sequenceNumber = itr->value[i]["i"].GetInt();
            tuple.voltage = itr->value[i]["voltage"].GetDouble();
            tuple.current = itr->value[i]["current"].GetDouble();
                inputBatch_Query2[(tuple.sequenceNumber / 1000) % queueSize_Query2].push_back(tuple);
        }

        for (int i = 0; i < queueSize_Query2; i++) {
            if (!inputBatch_Query2[i].empty()) {
                int batchSequence = inputBatch_Query2[i].back().sequenceNumber / 1000;
                if (inputBatch_Query2[i].size() == maxThrehold_Query2 || (batchNumber - batchSequence ) >= lateTime_Query2) {
                    double voltage[1000];
                    double current[1000];
                    for (int j = 0; j < 1000; j++) {
                        voltage[j] = 2;
                        current[j] = 2;
                    }
                    int tempQueueSize = inputBatch_Query2[i].size();
                    for (int j = 0; j < tempQueueSize; j++) {
                        Tuple tuple = inputBatch_Query2[i][j];
                        voltage[j] = tuple.voltage;
                        current[j] = tuple.current;
                    }
                    inputBatch_Query2[i].clear();
                    Feature X = incremental_Detector_Disorder.compute_input_signal(voltage, current, period_Query2,true);
                    X.sequenceNumber = batchSequence;
                    if(tempQueueSize >= 900) {
                        X.badPoint = -1;
                    }
                    else{
                        X.badPoint = 1;
                    }

                    vector<Result> resultBucket;
                    int eventPredict = incremental_Detector_Disorder.predictEvent(X,Radius_Query2);
                    vector<Event> event_interval_indices_vector = incremental_Detector_Disorder.predict(X.sequenceNumber, X);
                    if(eventPredict == 1){
                        waitBatch.push_back(X.sequenceNumber);
                    }
                    else{
                        Result result_temp;
                        result_temp.s = X.sequenceNumber;
                        result_temp.d = 0;
                        result_temp.event_s = -1;
                        resultBucket.push_back(result_temp);
                    }
                    waitBatch.sort();

                    if(!event_interval_indices_vector.empty()) {
                        for (int eventId = 0; eventId < event_interval_indices_vector.size(); eventId++) {
                            Event event_interval_indices = event_interval_indices_vector[eventId];
                            if (event_interval_indices.event_end != 0) {  // if an event is returned
                                cout << "batchNumber" <<  X.sequenceNumber << "  Event Detected at   " << event_interval_indices.event_start << "," <<
                                     event_interval_indices.event_end << endl;
                                int mean_event_index = (event_interval_indices.event_start + event_interval_indices.event_end) / 2 + 1;

                                Result result_temp;
                                result_temp.s = event_interval_indices.sequenceNumber;
                                result_temp.d = 1;
                                result_temp.event_s = mean_event_index;
                                resultBucket.push_back(result_temp);
                                waitBatch.remove(event_interval_indices.sequenceNumber);
                            }
                        }
                    }else{
                        incremental_Detector_Disorder.refresh_window(MAXIMUM_WINDOW_SIZE_Query2);
                    }
                    int tempWaitSize = waitBatch.size();
                    for(int i = 0; i< tempWaitSize; i++){
                        if(waitBatch.front() <= incremental_Detector_Disorder.mapFeatureBlockedPos){
                            Result result_temp;
                            result_temp.s = waitBatch.front();
                            result_temp.d = 0;
                            result_temp.event_s = -1;
                            resultBucket.push_back(result_temp);
                            waitBatch.pop_front();
                        }
                    }
                    event_interval_indices_vector.clear();

                    for(int i = 0; i< resultBucket.size(); i++){
                        std::string jsonTemp =   "{\"s\":";
                        jsonTemp.append(to_string (resultBucket[i].s));
                        if(resultBucket[i].d == 0){
                            jsonTemp.append(",\"d\":0,\"event_s\":");
                        }else{
                            jsonTemp.append(",\"d\":1,\"event_s\":");
                        }
                        if(resultBucket[i].event_s == -1){
                            jsonTemp.append(to_string(-1));
                        }else {
                            jsonTemp.append(to_string((int) resultBucket[i].event_s));
                        }
                        jsonTemp.append("}");


                        std::string response_string_post;
                        struct curl_slist *headers = NULL;
                        headers = curl_slist_append(headers, "Content-Type: application/json");
                        curl_easy_setopt(curl_post, CURLOPT_HTTPHEADER, headers);
                        curl_easy_setopt(curl_post, CURLOPT_POSTFIELDS, jsonTemp.c_str());
                        curl_easy_setopt(curl_post, CURLOPT_WRITEFUNCTION, writeFunction_Parallel2);
                        curl_easy_setopt(curl_post, CURLOPT_WRITEDATA, &response_string_post);
                        curl_easy_perform(curl_post);
                    }
                    resultBucket.clear();

                }
            }
        }
        batchNumber++;

    }
    cout << "> Query 2 done!" << endl;

}