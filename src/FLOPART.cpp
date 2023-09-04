#include <vector>
#include <stdio.h>
#include "funPieceListLog.h"
#include <math.h>
#include <R.h> // for Rprintf
#include "FLOPART.h"

class Rule {
public:
  void update
  (PiecewisePoissonLossLog *this_cost,
   PiecewisePoissonLossLog *this_cost_prev,
   PiecewisePoissonLossLog *other_cost_prev,
   PiecewisePoissonLossLog *initial_cost,
   int curr_label_type,
   bool at_beginning,
   bool at_end,
   double cum_weight_prev_i,
   double cum_weight_i,
   const double *weight_vec,
   int data_i,
   const int *data_vec,
   int verbose,
   double penalty
   ){
    PiecewisePoissonLossLog cost_of_change;
    if(data_i == 0){
      *this_cost = *initial_cost;
    } else if(no_change(curr_label_type, at_beginning, at_end)){
      *this_cost = *this_cost_prev;
    } else if(infinite(curr_label_type, at_beginning, at_end)){
      this_cost->set_infinite();
    } else if(change(curr_label_type, at_beginning, at_end)){
      min_more_or_less(&cost_of_change, other_cost_prev, verbose);
      cost_of_change.set_prev_seg_end(data_i-1);
      cost_of_change.addPenalty(penalty, cum_weight_prev_i);
      this_cost->set_to_min_env_of(this_cost_prev, &cost_of_change, verbose);
    }
    this_cost->addDataLoss
      (cum_weight_prev_i, cum_weight_i, weight_vec,
       data_i, data_vec);
  }
  virtual void min_more_or_less
  (PiecewisePoissonLossLog*, PiecewisePoissonLossLog*, int) = 0;
  virtual bool no_change(int,bool,bool) = 0;
  virtual bool infinite(int,bool,bool) = 0;
  virtual bool change(int,bool,bool) = 0;
};

class Up : public Rule {
  void min_more_or_less
  (PiecewisePoissonLossLog *cost_of_change_up,
   PiecewisePoissonLossLog *down_cost_prev,
   int verbose){
    cost_of_change_up->set_to_min_less_of(down_cost_prev, verbose);
  }
  bool no_change(int curr_label_type, bool at_beginning, bool at_end){
    return curr_label_type == LABEL_PEAKEND && !at_beginning && !at_end;
  }
  bool infinite(int curr_label_type, bool at_beginning, bool at_end){
    return curr_label_type == LABEL_NOPEAKS
      || (curr_label_type == LABEL_PEAKSTART && at_beginning)
      || (curr_label_type == LABEL_PEAKEND && at_end);
  }
  bool change(int curr_label_type, bool at_beginning, bool at_end){
    return curr_label_type == LABEL_UNLABELED 
      || (curr_label_type== LABEL_PEAKSTART && !at_beginning) 
      || (curr_label_type == LABEL_PEAKEND && at_beginning);
  }
};
    
class Down : public Rule {
  void min_more_or_less
  (PiecewisePoissonLossLog *cost_of_change_down,
   PiecewisePoissonLossLog *up_cost_prev,
   int verbose){
    cost_of_change_down->set_to_min_more_of(up_cost_prev, verbose);
  }
  bool no_change(int curr_label_type, bool at_beginning, bool at_end){
    return (curr_label_type == LABEL_NOPEAKS && !at_beginning)
      || (curr_label_type == LABEL_PEAKSTART && !at_beginning && !at_end);
  }
  bool infinite(int curr_label_type, bool at_beginning, bool at_end){
    return (curr_label_type == LABEL_PEAKEND && at_beginning) ||
      (curr_label_type == LABEL_PEAKSTART && at_end);
  }
  bool change(int curr_label_type, bool at_beginning, bool at_end){
    return curr_label_type == LABEL_UNLABELED 
      || (curr_label_type == LABEL_PEAKSTART && at_beginning)
      || (curr_label_type == LABEL_NOPEAKS && at_beginning)
      || (curr_label_type == LABEL_PEAKEND && !at_beginning);
  }
};

int FLOPART
  (const int *data_vec, const double *weight_vec, const int data_count,
   const double penalty,
   const int *label_types,
   const int *label_starts,
   const int *label_ends,
   const int label_count,
   // the following matrices are for output.
   // cost_mat and intervals_mat store the optimal cost and number of intervals
   // at each time point, for the up and down cost models.
   // end_vec and mean_vec store the best model up to and including the
   // last data point. 
   double *cost_mat, //data_count x 2.
   int *end_vec,//data_count
   double *mean_vec,//data_count
   int *intervals_mat,//data_count x 2
   int *state_vec
   ){
  for(int label_i=0; label_i < label_count; label_i++){
    if(label_starts[label_i] < 0){
      return ERROR_LABEL_START_MUST_BE_AT_LEAST_ZERO;
    }
    if(data_count <= label_ends[label_i]){
      return ERROR_LABEL_END_MUST_BE_LESS_THAN_DATA_SIZE;
    }
    if(label_ends[label_i] <= label_starts[label_i]){
      return ERROR_LABEL_END_MUST_BE_GREATER_THAN_LABEL_START;
    }
    if(0 < label_i && label_starts[label_i] <= label_ends[label_i-1]){
      return ERROR_LABEL_START_SHOULD_BE_GREATER_THAN_PREVIOUS_LABEL_END;
    }
    switch(label_types[label_i]){
    case LABEL_NOPEAKS:
    case LABEL_PEAKSTART:
    case LABEL_PEAKEND:
      break;
    default:
      Rprintf("label[%d]=%d\n", label_i, label_types[label_i]);
      return ERROR_UNRECOGNIZED_LABEL_TYPE;
    }
  }
  bool at_beginning = false;
  bool at_end = false;  
  int curr_label_index = 0;
  int curr_label_type = LABEL_UNLABELED; //1 = peakStart, 0 = noPeak, -1 = peakEnd
  double min_log_mean=INFINITY, max_log_mean=-INFINITY;
  for(int data_i=0; data_i<data_count; data_i++){
    
    double log_data = log( (double)(data_vec[data_i]) );
    if(log_data < min_log_mean){
      min_log_mean = log_data;
    }
    if(max_log_mean < log_data){
      max_log_mean = log_data;
    }
  }
  if(min_log_mean == max_log_mean){
    return ERROR_MIN_MAX_SAME;
  }
  CostMatrix cost_funs(data_count);
  PiecewisePoissonLossLog
    *up_cost, *down_cost, *up_cost_prev, *down_cost_prev;
  // initialize to avoid warnings.
  up_cost_prev = &cost_funs.cost_vec[0];
  down_cost_prev = &cost_funs.cost_vec[0];
  PiecewisePoissonLossLog initial_cost;
  // initialization C_1(m)=loss_1(m)/w_1
  initial_cost.piece_list.emplace_back
    (1.0, -data_vec[0], 0.0,
     min_log_mean, max_log_mean, -1, false);
  int verbose=0;
  double cum_weight_i = 0.0, cum_weight_prev_i = 0.0;
  for(int data_i=0; data_i<data_count; data_i++){
    cum_weight_i += weight_vec[data_i];
    up_cost = &cost_funs.cost_vec[data_i];
    down_cost = &cost_funs.cost_vec[data_i + data_count];
    at_beginning = false;
    at_end = false;
    //if we haven't checked all the labels 
    if(curr_label_index < label_count){
      if(data_i >= label_starts[curr_label_index] &&
	 data_i <= label_ends[curr_label_index]){
        //if >= current start and <= current end
        curr_label_type = label_types[curr_label_index];
        if(data_i == label_starts[curr_label_index]){
          at_beginning = true;
        }
        if(data_i == label_ends[curr_label_index]){
          at_end = true;
          curr_label_index++;
        }
      }else{
        curr_label_type = LABEL_UNLABELED;
      }
    }else{
      curr_label_type = LABEL_UNLABELED;
    }
    for(int state=0; state<=1; state++){
      if(state==0){
	Down down;
	down.update
	  (down_cost, down_cost_prev, up_cost_prev,
	   &initial_cost, curr_label_type, at_beginning, at_end,
	   cum_weight_prev_i, cum_weight_i, weight_vec,
	   data_i, data_vec, verbose, penalty);
      }else{
	Up up;
	up.update
	  (up_cost, up_cost_prev, down_cost_prev,
	   &initial_cost, curr_label_type, at_beginning, at_end,
	   cum_weight_prev_i, cum_weight_i, weight_vec,
	   data_i, data_vec, verbose, penalty);
      }
    }
    cum_weight_prev_i = cum_weight_i;
    up_cost_prev = up_cost;
    down_cost_prev = down_cost;
  }
  // Then store number of intervals and optimal cost value, for each
  // cost function.
  cost_funs.copy_min_cost_intervals(cost_mat, intervals_mat);
  double cost = cost_funs.decode_optimal_mean_end_state
    (mean_vec, end_vec, state_vec);
  if(cost == INFINITY)return ERROR_INFINITE_COST;
  return 0;
}

