#include <Rcpp.h>
#include <R.h>
#include "FLOPART.h"

//' Lookup the integer values used to represent different label types
//'
//' @return Integer vector with names corresponding to supported label types
// [[Rcpp::export]]
Rcpp::IntegerVector get_label_code(){
  Rcpp::IntegerVector code {LABEL_NOPEAKS, LABEL_PEAKSTART, LABEL_PEAKEND};
  code.names() = Rcpp::CharacterVector {"noPeaks", "peakStart", "peakEnd"};
  return code;
}

//' Interface to FLOPART C++ code
//'
//' @param data_vec Integer vector of non-negative count data
//' @param weight_vec Numeric vector of positive weights (same size as data_vec)
//' @param penalty non-negative real-valued penalty (larger for fewer peaks)
//' @param label_type_vec Integer vector of label types
//' @param label_start_vec Integer vector of label starts
//' @param label_end_vec Integer vector of label ends
//' @return List with named elements: cost_mat and intervals_mat (one row for each data point, first column up, second down), segments_df (one row for each segment in the optimal model)
// [[Rcpp::export]]
Rcpp::List FLOPART_interface
(const Rcpp::IntegerVector data_vec,
 const Rcpp::NumericVector weight_vec,
 const double penalty,
 const Rcpp::IntegerVector label_type_vec,
 const Rcpp::IntegerVector label_start_vec,
 const Rcpp::IntegerVector label_end_vec){
  int data_count = data_vec.size();
  if(data_count < 2){
    Rcpp::stop("need at least two data points");
  }
  if(data_count != weight_vec.size()){
    Rcpp::stop("data_vec and weight_vec should be same size");
  }
  if(Rcpp::any(!Rcpp::is_finite(data_vec))){
    Rcpp::stop("data must be finite");
  }
  int label_count = label_type_vec.size();
  if(label_count != label_start_vec.size()){
    Rcpp::stop("label_start_vec and label_type_vec should be same size");
  }
  if(label_count != label_end_vec.size()){
    Rcpp::stop("label_end_vec and label_type_vec should be same size");
  }
  Rcpp::NumericMatrix cost_mat(data_count, 2);
  Rcpp::IntegerMatrix intervals_mat(data_count, 2);
  Rcpp::IntegerVector rev_end_vec(data_count);
  Rcpp::NumericVector rev_mean_vec(data_count);
  Rcpp::IntegerVector rev_state_vec(data_count);
  const int *label_type_ptr=0, *label_start_ptr=0, *label_end_ptr=0;
  if(0 < label_count){
    label_type_ptr = &label_type_vec[0];
    label_start_ptr = &label_start_vec[0];
    label_end_ptr = &label_end_vec[0];
  }
  int status = FLOPART
    (&data_vec[0], &weight_vec[0],
     data_count, penalty,
     label_type_ptr, label_start_ptr, label_end_ptr, label_count,
     &cost_mat[0], &rev_end_vec[0], &rev_mean_vec[0], &intervals_mat[0],
     &rev_state_vec[0]
     );
  if(status == ERROR_MIN_MAX_SAME){
    Rcpp::stop("data[i]=%d for all i, fix by using a data sequence with more than one value", data_vec[0]);
  }
  if(status == ERROR_LABEL_END_MUST_BE_GREATER_THAN_LABEL_START){
    Rcpp::stop("label end must be greater than label start, use FLOPART_data(counts, labels) to see firstRow/lastRow values passed to C++ code, and fix by removing labels with firstRow==lastRow");
  }
  if(status == ERROR_LABEL_START_SHOULD_BE_GREATER_THAN_PREVIOUS_LABEL_END){
    Rcpp::stop("label start should be greater than previous label end, use FLOPART_data(counts, labels) to see firstRow/lastRow columns passed to C++ code");
  }
  if(status == ERROR_UNRECOGNIZED_LABEL_TYPE){
    Rcpp::stop("unrecognized label type, use FLOPART_data(counts, labels) to see type column passed to C++ code, and get_label_code() for a mapping between annotation/character and type/integer values");
  }
  if(status == ERROR_LABEL_END_MUST_BE_LESS_THAN_DATA_SIZE){
    Rcpp::stop("label end must be less than data size, use FLOPART_data(counts, labels) to see lastRow column passed to C++ code");
  }
  if(status == ERROR_LABEL_START_MUST_BE_AT_LEAST_ZERO){
    Rcpp::stop("label start must be at least zero, use FLOPART_data(counts, labels) to see firstRow column passed to C++ code");
  }
  Rcpp::DataFrame seg_df;
  if(status == ERROR_INFINITE_COST){
    seg_df = Rcpp::DataFrame::create();
  }else{
    //convert to segments data table.
    int seg_count=1;
    while(seg_count < data_count && 0 <= rev_end_vec[seg_count-1]){
      seg_count++;
    }
    Rcpp::NumericVector seg_mean_vec(seg_count);
    Rcpp::IntegerVector seg_start_vec(seg_count);
    Rcpp::IntegerVector seg_end_vec(seg_count);
    Rcpp::IntegerVector seg_state_vec(seg_count);
    for(int seg_i=0; seg_i < seg_count; seg_i++){
      int mean_index = seg_count-1-seg_i;
      seg_mean_vec[seg_i] = rev_mean_vec[mean_index];
      seg_state_vec[seg_i] = rev_state_vec[mean_index];
      if(mean_index==0){
        seg_end_vec[seg_i] = data_count;
      }else{
        seg_end_vec[seg_i] = rev_end_vec[mean_index-1]+1;
      }
      if(seg_i==0){
        seg_start_vec[seg_i] = 1;
      }else{
        seg_start_vec[seg_i] = rev_end_vec[mean_index]+2;
      }
    }
    seg_df = Rcpp::DataFrame::create
      (Rcpp::Named("mean", seg_mean_vec),
       Rcpp::Named("firstRow", seg_start_vec),
       Rcpp::Named("lastRow", seg_end_vec),
       Rcpp::Named("state", seg_state_vec)
       );
  }
  return Rcpp::List::create
    (Rcpp::Named("cost_mat", cost_mat),
     Rcpp::Named("intervals_mat", intervals_mat),
     Rcpp::Named("segments_df", seg_df)
     );
}
