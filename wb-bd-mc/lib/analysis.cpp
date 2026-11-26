#include "analysis.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {
double correlate(MATRIX& input,long column,long threads,MATRIX& result,long output_column,
                 POST_PROCESSING::CorrelationWindow window,bool center) {
    using Window=POST_PROCESSING::CorrelationWindow;
    if(threads<=0 || input.rows<=0 || input.cols<=0 || column<0 || column>=input.cols ||
       result.rows<=0 || output_column<0 || output_column>=result.cols || &input==&result)
        throw std::invalid_argument("invalid autocorrelation shape, column, alias or thread count");
    if(window!=Window::LegacyFixedWindow && window!=Window::AllOrigins)
        throw std::invalid_argument("unknown autocorrelation window");
    // H-D11: fixed-window and all-origin estimators are different valid choices.
    // Normalize by the chosen origins at each lag; do not silently reinterpret
    // a short window as all origins or include circular wraparound pairs.
    const long lags=window==Window::LegacyFixedWindow?input.rows/2:result.rows;
    if(lags<=0 || result.rows<lags || lags>input.rows ||
       (window==Window::LegacyFixedWindow && input.rows<4))
        throw std::invalid_argument("autocorrelation has no valid origins or insufficient result rows");
    long double mean=0.;
    for(long i=0;i<input.rows;++i) {
        const double value=input(i,column);
        if(!std::isfinite(value)) throw std::invalid_argument("autocorrelation input must be finite");
        if(center) mean+=static_cast<long double>(value)/input.rows;
    }
    std::vector<double> values(static_cast<std::size_t>(lags));
#pragma omp parallel for default(none) if(threads > 1) num_threads(threads) \
    shared(input, column, lags, window, mean, values)
    for(long lag=0;lag<lags;++lag) {
        const long origins=window==Window::AllOrigins?input.rows-lag:input.rows/2-1;
        long double sum=0.;
        for(long origin=0;origin<origins;++origin)
            sum+=(input(origin,column)-mean)*(input(origin+lag,column)-mean);
        values[static_cast<std::size_t>(lag)]=static_cast<double>(sum/origins);
    }
    for(double v:values) if(!std::isfinite(v)) throw std::overflow_error("autocorrelation result overflows");
    for(long lag=0;lag<lags;++lag) result(lag,output_column)=values[static_cast<std::size_t>(lag)];
    return 0.;
}
std::vector<std::string> tokens(std::string line) {
    for(char& c:line) if(c==',' || c=='\t') c=' ';
    std::istringstream stream(line); std::vector<std::string> out;
    for(std::string token;stream>>token;) out.push_back(token);
    return out;
}
long nonnegative(const std::string& value,const char* what) {
    std::size_t used=0;
    long out=0;
    try {out=std::stol(value,&used);} catch(const std::exception&) {throw std::invalid_argument(what);}
    if(used!=value.size() || out<0) throw std::invalid_argument(what);
    return out;
}
bool numeric(const std::string& value,double& out) {
    try {std::size_t used=0;out=std::stod(value,&used);return used==value.size();}
    catch(const std::exception&) {return false;}
}
long resolve(const std::string& name,const std::vector<std::string>& header,long columns) {
    long index=-1;
    if(!name.empty() && std::all_of(name.begin(),name.end(),[](char c){return c>='0' && c<='9';}))
        index=nonnegative(name,"column index is out of range");
    else {
        const auto item=std::find(header.begin(),header.end(),name);
        if(item==header.end() || header.size()!=static_cast<std::size_t>(columns))
            throw std::invalid_argument("named column requires a matching table header: "+name);
        index=static_cast<long>(item-header.begin());
    }
    if(index<0 || index>=columns) throw std::invalid_argument("column is out of range: "+name);
    return index;
}
void help(bool connector) {
    std::cout<<"Usage: "<<(connector?"correlation_function_individual_chain":"measure_stress_autocorrelation")
        <<" INPUT OUTPUT --time-column NAME_OR_ZERO_BASED_INDEX --columns C1,C2,...\n"
          "  [--window all_origins|legacy_fixed_window] [--max-lag NONNEGATIVE]\n"
          "  [--center] [--threads POSITIVE]\n"
          "Whitespace/CSV table; optional header or # header; uniform increasing time.\n";
    if(connector) std::cout<<"Connector mode: explicit component columns are triples (qx,qy,qz) per chain.\n"
        "Output is mean chain dot-product correlation, followed by each chain.\n"
        "This reconstructed tool does not interpret legacy .chain files.\n";
}
} // namespace

double POST_PROCESSING::compute_autocorrelation_brute_force(MATRIX& input,const long column,
    MATRIX& result,const long output_column,CorrelationWindow window,bool center) {
    return correlate(input,column,1,result,output_column,window,center);
}
double POST_PROCESSING::compute_autocorrelation_brute_force_OMP(MATRIX& input,const long column,
    const long threads,MATRIX& result,const long output_column,CorrelationWindow window,bool center) {
    return correlate(input,column,threads,result,output_column,window,center);
}

int POST_PROCESSING::autocorrelation_cli(int argc,char* argv[],bool connector) {
    if(argc==1 || (argc==2 && std::string(argv[1])=="--help")) {help(connector);return 0;}
    try {
        if(argc<3) throw std::invalid_argument("INPUT and OUTPUT are required");
        std::string time_name,window_name="all_origins";
        std::vector<std::string> component_names;
        long threads=1,max_lag=-1; bool center=false;
        for(int i=3;i<argc;++i) {
            const std::string flag=argv[i];
            if(flag=="--center") {center=true;continue;}
            if(i+1>=argc) throw std::invalid_argument("option requires a value: "+flag);
            const std::string value=argv[++i];
            if(flag=="--time-column") time_name=value;
            else if(flag=="--columns") component_names=tokens(value);
            else if(flag=="--threads") threads=nonnegative(value,"threads must be a positive integer");
            else if(flag=="--max-lag") max_lag=nonnegative(value,"max-lag must be a nonnegative integer");
            else if(flag=="--window") window_name=value;
            else throw std::invalid_argument("unknown option: "+flag);
        }
        if(time_name.empty() || component_names.empty() || threads<=0)
            throw std::invalid_argument("explicit time/component columns and positive threads are required");
        if(connector && component_names.size()%3!=0)
            throw std::invalid_argument("connector component columns must form complete qx,qy,qz triples");
        if(window_name!="all_origins" && window_name!="legacy_fixed_window")
            throw std::invalid_argument("window must be all_origins or legacy_fixed_window");
        const auto window=window_name=="all_origins"?CorrelationWindow::AllOrigins:CorrelationWindow::LegacyFixedWindow;
        if(std::filesystem::exists(argv[2])) throw std::invalid_argument("output already exists; choose a fresh path");
        std::ifstream input(argv[1]);
        if(!input) throw std::invalid_argument("cannot open input table");
        std::vector<double> raw; std::vector<std::string> header;
        long rows=0,columns=0,line_number=0;
        for(std::string line;std::getline(input,line);) {
            ++line_number;
            const auto first=line.find_first_not_of(" \t\r");
            if(first==std::string::npos) continue;
            if(line[first]=='#') {
                if(rows==0) {
                    const auto candidate=tokens(line.substr(first+1));
                    if(!candidate.empty()) header=candidate;
                }
                continue;
            }
            const auto fields=tokens(line); if(fields.empty()) continue;
            double first_value=0.;
            if(rows==0 && !numeric(fields.front(),first_value)) {header=fields;continue;}
            if(columns==0) columns=static_cast<long>(fields.size());
            if(fields.size()!=static_cast<std::size_t>(columns))
                throw std::invalid_argument("ragged input at line "+std::to_string(line_number));
            for(const auto& field:fields) {
                double value=0.;
                if(!numeric(field,value) || !std::isfinite(value))
                    throw std::invalid_argument("non-finite/nonnumeric input at line "+std::to_string(line_number));
                raw.push_back(value);
            }
            ++rows;
        }
        if(rows<2) throw std::invalid_argument("at least two samples are required");
        const long time_column=resolve(time_name,header,columns);
        std::vector<long> selected;std::set<long> unique;
        for(const auto& name:component_names) {
            const long index=resolve(name,header,columns);
            if(index==time_column || !unique.insert(index).second)
                throw std::invalid_argument("component columns must be distinct and exclude time");
            selected.push_back(index);
        }
        const double dt=raw[static_cast<std::size_t>(columns+time_column)]-raw[static_cast<std::size_t>(time_column)];
        if(!std::isfinite(dt) || dt<=0) throw std::invalid_argument("time must be strictly increasing and uniform");
        for(long r=1;r<rows;++r) {
            const double delta=raw[static_cast<std::size_t>(r*columns+time_column)]-raw[static_cast<std::size_t>((r-1)*columns+time_column)];
            const double tolerance=1e-9*std::max(std::abs(dt),std::abs(delta))+
                32*std::numeric_limits<double>::epsilon()*std::abs(raw[static_cast<std::size_t>(r*columns+time_column)]);
            if(delta<=0 || !std::isfinite(delta) || std::abs(delta-dt)>tolerance)
                throw std::invalid_argument("time must be strictly increasing and uniform");
        }
        if(max_lag<0) max_lag=rows/2-1;
        if(max_lag>=rows || (window==CorrelationWindow::LegacyFixedWindow && (rows<4 || max_lag>=rows/2)))
            throw std::invalid_argument("max-lag is outside the selected estimator's valid range");
        const long computed_rows=window==CorrelationWindow::LegacyFixedWindow?rows/2:max_lag+1;
        MATRIX series(rows,static_cast<long>(selected.size()),0.);
        for(long r=0;r<rows;++r) for(std::size_t c=0;c<selected.size();++c)
            series(r,static_cast<long>(c))=raw[static_cast<std::size_t>(r*columns+selected[c])];
        MATRIX correlations(computed_rows,static_cast<long>(selected.size()),0.);
        for(std::size_t c=0;c<selected.size();++c)
            compute_autocorrelation_brute_force_OMP(series,static_cast<long>(c),threads,correlations,static_cast<long>(c),window,center);
        std::ofstream output(argv[2]);
        if(!output) throw std::runtime_error("cannot open output table");
        output<<"# estimator="<<window_name<<" centered="<<(center?"true":"false")<<" zero_based_columns=true\n";
        output<<"# lag_time origin_count";
        if(connector) {
            output<<" mean_chain_dot";
            for(std::size_t c=0;c<selected.size()/3;++c) output<<" chain_"<<c<<"_dot";
        } else for(const auto& name:component_names) output<<" C_"<<name;
        output<<'\n'<<std::setprecision(17);
        for(long lag=0;lag<=max_lag;++lag) {
            output<<lag*dt<<' '<<(window==CorrelationWindow::AllOrigins?rows-lag:rows/2-1);
            if(connector) {
                std::vector<long double> dot(selected.size()/3,0.);
                long double mean=0.;
                for(std::size_t c=0;c<dot.size();++c) {
                    for(long k=0;k<3;++k) dot[c]+=correlations(lag,static_cast<long>(3*c)+k);
                    mean+=dot[c]/dot.size();
                }
                output<<' '<<mean;for(auto value:dot) output<<' '<<value;
            } else for(std::size_t c=0;c<selected.size();++c) output<<' '<<correlations(lag,static_cast<long>(c));
            output<<'\n';
        }
        if(!output) throw std::runtime_error("failed writing correlation output");
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<"autocorrelation: "<<error.what()<<'\n';return 2;
    }
}
