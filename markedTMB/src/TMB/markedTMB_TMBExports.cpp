// Generated by TMBtools: do not edit by hand

#define TMB_LIB_INIT R_init_markedTMB_TMBExports
#include <TMB.hpp>
#include "msjs.hpp"
#include "msjsu.hpp"
#include "msld.hpp"
#include "multistate.hpp"
#include "mvms.hpp"
#include "smsld.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if(model == "msjs") {
    return msjs(this);
  } else if(model == "msjsu") {
    return msjsu(this);
  } else if(model == "msld") {
    return msld(this);
  } else if(model == "multistate") {
    return multistate(this);
  } else if(model == "mvms") {
    return mvms(this);
  } else if(model == "smsld") {
    return smsld(this);
  } else {
    error("Unknown model.");
  }
  return 0;
}