/*
  ____  _            _           _
 | __ )| | ___   ___| | __   ___| | __ _ ___ ___
 |  _ \| |/ _ \ / __| |/ /  / __| |/ _` / __/ __|
 | |_) | | (_) | (__|   <  | (__| | (_| \__ \__ \
 |____/|_|\___/ \___|_|\_\  \___|_|\__,_|___/___/

*/
#include "block_la.h"
#include "defines.h"
#include <ctype.h>     // toupper
#include <math.h>      // pow
#include <sys/param.h> // MIN

/*
  ____        __ _       _ _   _
 |  _ \  ___ / _(_)_ __ (_) |_(_) ___  _ __  ___
 | | | |/ _ \ |_| | '_ \| | __| |/ _ \| '_ \/ __|
 | |_| |  __/  _| | | | | | |_| | (_) | | | \__ \
 |____/ \___|_| |_|_| |_|_|\__|_|\___/|_| |_|___/

*/

// Velocity profile data
typedef struct {
  data_t a, d;             // actual accelerations
  data_t f, l;             // actual feedrate and length
  data_t fs, fe;           // initial and final feedrate
  data_t dt_1, dt_m, dt_2; // durations
  data_t dt;               // total duration

/***************************************************************************************************/
// LOOK_AHEAD Implementation
  data_t v_i, v_m, v_f;
  data_t ds_1, ds_m, ds_2;  // abscissa
/***************************************************************************************************/
} block_profile_t;

typedef struct block {
  char *line;               // G-code line as a string
  block_type_t type;        // block type
  size_t n;                 // block number
  size_t tool;              // tool number
  data_t feedrate;          // nominal feedrate
  data_t arc_feedrate;      // actual feedrate along an arc
  data_t spindle;           // spindle rotational speed
  point_t *target;          // target position
  point_t *delta;           // segment projections
  point_t *center;          // arc center coordinates
  data_t length;            // segment/arc length
  data_t i, j, r;           // arc parameters (offsets and radius)
  data_t theta_0, dtheta;   // initial angle and arc angle
  data_t acc;               // actual acceleration
  machine_t const *machine; // the machine reference
  block_profile_t *prof;    // the speed profile
  struct block *next;       // reference to the next block
  struct block *prev;       // reference to the previous block
} block_t;

/*
  _____                 _   _
 |  ___|   _ _ __   ___| |_(_) ___  _ __  ___
 | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
 |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
 |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/

*/

/* STATIC FUNCTIONS ***********************************************************/
static point_t *start_point(block_t const *b);
static ccnc_error_t block_set_fields(block_t *b, char cmd, char *arg);
// static void block_compute(block_t *b);
static ccnc_error_t block_arc(block_t *b);
static data_t quantize(data_t t, data_t tq, data_t *dq);
static ccnc_error_t block_parse(block_t *b);

/***************************************************************************************************/
// LOOK_AHEAD Implementation
static data_t curvlin_to_time(data_t v_0, data_t a, data_t ds);
/***************************************************************************************************/

/* LIFECYCLE ******************************************************************/
block_t *block_new(char const *line, block_t *prev, machine_t const *machine) {
  assert(line);
  block_t *b = (block_t *)malloc(sizeof(block_t));
  if (!b) {
    eprintf("Could not allocate memory for block line %s\n", line);
    goto fail;
  }

  if (prev) { // this is not the first block, copy memory from previous one
    memcpy(b, prev, sizeof(*b));
    b->prev = prev;
    b->next = NULL;
    prev->next = b;
  } else { // this is the very first block
    memset(b, 0, sizeof(*b));
  }

  b->machine = machine;
  b->acc = machine_A(machine);

  b->prof = (block_profile_t *)malloc(sizeof(block_profile_t));
  if (!b->prof) {
    eprintf("Could not allocate memory for velocity profile\n");
    goto fail;
  }
  b->target = point_new();
  b->delta = point_new();
  b->center = point_new();
  if (!b->center || !b->delta || !b->target) {
    eprintf("Could not allocate memory for block points\n");
    goto fail;
  }
  b->line = strdup(line);
  if (!b->line) {
    eprintf("Could not allcate memory for G-code string\n");
    goto fail;
  }

  // Clear i, j & r for double set flag
  if (b->i) b->i = 0;
  if (b->j) b->j = 0;
  if (b->r) b->r = 0;

  if (block_parse(b) != NO_ERR) {
    eprintf("Could not parse block\n");
    goto fail;
  }

  return b;

fail:
  if (b)
    block_free(b);
  return NULL;
}

void block_free(block_t *b) {
  assert(b);
  if (b->line)
    free(b->line);
  if (b->prof)
    free(b->prof);
  if (b->target)
    point_free(b->target);
  if (b->center)
    point_free(b->center);
  if (b->delta)
    point_free(b->delta);
  free(b);
  b = NULL;
}

void block_print(block_t const *b, FILE *out) {
  assert(b && out);
  char *start = NULL, *end = NULL;
  point_t *p0 = start_point(b);
  point_inspect(p0, &start);
  point_inspect(b->target, &end);
  fprintf(out, "%03lu %s->%s F%7.1f S%7.1f T%02lu G%02d\n", b->n, start, end,
          b->feedrate, b->spindle, b->tool, b->type);
  free(start);
  free(end);
}

/* ACCESSORS ******************************************************************/

// #define block_getter(typ, par, name)                                           \
//   typ block_##name(block_t const *b) {                                         \
//     assert(b);                                                                 \
//     return b->par;                                                             \
//   }

// block_getter(data_t, length, length);
// block_getter(data_t, dtheta, dtheta);
// block_getter(block_type_t, type, type);
// block_getter(data_t, prof->dt, dt);
// block_getter(char *, line, line);
// block_getter(size_t, n, n);
// block_getter(data_t, r, r);
// block_getter(point_t *, center, center);
// block_getter(point_t *, target, target);
// block_getter(block_t *, next, next);

class_getter(block, block_t, data_t, length, length);
class_getter(block, block_t, data_t, dtheta, dtheta);
class_getter(block, block_t, block_type_t, type, type);
class_getter(block, block_t, data_t, prof->dt, dt);
class_getter(block, block_t, char *, line, line);
class_getter(block, block_t, size_t, n, n);
class_getter(block, block_t, data_t, r, r);
class_getter(block, block_t, point_t *, center, center);
class_getter(block, block_t, point_t *, target, target);
class_getter(block, block_t, block_t *, next, next);
class_getter(block, block_t, block_t *, prev, prev);

/* METHODS ********************************************************************/

/***************************************************************************************************/
// LOOK_AHEAD Implementation

// void block_velocities(block_t *b, data_t t, data_t *s) {
//   assert(b);
//   data_t r;
//   data_t dt_1 = b->prof->dt_1;
//   data_t dt_2 = b->prof->dt_2;
//   data_t dt_m = b->prof->dt_m;
//   data_t a = b->prof->a;
//   data_t d = b->prof->d;
//   data_t f = b->prof->f;

//   // if (t < 0) {
//   //   exit; 
//   // }
//   //  else if 
// }
/***************************************************************************************************/

data_t block_lambda(block_t *b, data_t t, data_t *s) {
  // s - speed (feedrate)
  assert(b);
  data_t r;
  data_t dt_1 = b->prof->dt_1;
  data_t dt_2 = b->prof->dt_2;
  data_t dt_m = b->prof->dt_m;
  data_t a = b->prof->a;
  data_t d = b->prof->d;
  data_t f = b->prof->f;

  if (t < 0) {
    r = 0.0;
    *s = 0.0;
  } else if (t < dt_1) { // acceleration
    r = a * pow(t, 2) / 2.0;
    *s = a * t;
  } else if (t < dt_1 + dt_m) { // maintenance
    r = f * (dt_1 / 2.0 + (t - dt_1));
    *s = f;
  } else if (t < dt_1 + dt_m + dt_2) { // deceleration
    data_t t_2 = dt_1 + dt_m;
    r = f * dt_1 / 2.0 + f * (dt_m + t - t_2) +
        d / 2.0 * (pow(t, 2) + pow(t_2, 2)) - d * t * t_2;
    *s = f + d * (t - t_2);
  } else {
    r = b->prof->l;
    *s = 0.0;
  }
  
  r /= b->prof->l;
  *s *= 60; // convert to mm/min
  return r;
}

// /***************************************************************************************************/
// // LOOK_AHEAD Implementation
// point_t *block_interpolate(block_t *b, data_t lambda) 
// /***************************************************************************************************/

point_t *block_interpolate(block_t *b, data_t lambda) {
  assert(b);
  point_t *result = machine_setpoint(b->machine);
  point_t *p0 = start_point(b);
  
  // 1. the block describes a segment
  // x(t) = x(0) + d_x * lambda(t)
  // y(t) = y(0) + d_y * lambda(t)
  if (b->type == LINE) {
    point_set_x(result, point_x(p0) + point_x(b->delta) * lambda);
    point_set_y(result, point_y(p0) + point_y(b->delta) * lambda);
  }

  // 2. the block describes an arc
  // x(t) = x_c + R cos(theta_0 + dtheta * lambda(t))
  // y(t) = y_c + R sin(theta_0 + dtheta * lambda(t))
  else if (b->type == CWA || b->type == CCWA) {
    data_t angle = b->theta_0 + b->dtheta * lambda;
    point_set_x(result, point_x(b->center) + b->r * cos(angle));
    point_set_y(result, point_y(b->center) + b->r * sin(angle));
  } else {
    wprintf("Unexpected block type in interpolation\n");
    return NULL;
  }
  point_set_z(result, point_z(p0) + point_z(b->delta) * lambda);
  return result;
}

point_t *block_interpolate_t(block_t *b, data_t t, data_t *lambda, data_t *s) {
  *lambda = block_lambda(b, t, s);
  return block_interpolate(b, *lambda);
}

/* STATIC FUNCTIONS
 * ***********************************************************/

static point_t *start_point(block_t const *b) {
  assert(b);
  return b->prev ? b->prev->target : machine_zero(b->machine);
}

// "N01 G00 Z1000 Y500.10 T25 S5000 X123.321"
static ccnc_error_t block_parse(block_t *b) {
  assert(b);
  char *word, *line, *to_free;
  ccnc_error_t error = NO_ERR;
  point_t *p0;

  to_free = line = strdup(b->line);
  if (!line) {
    eprintf("Could not duplicate G-code line %s\n", b->line);
    return ALLOC_ERR;
  }

  // Tokenization
  // We shall remove any spaces at the end of the line before tokeinzing it
  while ((word = strsep(&line, " ")) != NULL) {
    // A valid command has at least two letters (e.g. G0)
    // This also skips double spaces or spaces at the end of the line
    if (strlen(word) < 2)
      continue;
    // word[0] is the G-code command
    // word + 1 is the argument
    error = block_set_fields(b, toupper(word[0]), word + 1);
    if (error != NO_ERR)
      break;
  }
  free(to_free);

  // Inherit coords from prev block
  p0 = start_point(b);
  point_modal(p0, b->target);
  point_delta(p0, b->target, b->delta);
  b->length = point_dist(p0, b->target);

  // Deal with motion blocks
  switch (b->type) {
  case LINE: // G01
    b->acc = machine_A(b->machine);
    b->arc_feedrate = b->feedrate;
    // block_compute(b);
    break;
  case CWA:  // G02
  case CCWA: // G03
    if (block_arc(b)) {
      wprintf("Could not calculate arc parameters\n");
      error = ARC_ERR;
      break;
    }
    b->arc_feedrate = MIN(
        b->feedrate,
        pow(3.0 / 4.0 * pow(machine_A(b->machine), 2) * pow(b->r, 2), 0.25) *
            60);
    // block_compute(b);
    break;
  default:
    break;
  }
  return error;
}

static ccnc_error_t block_set_fields(block_t *b, char cmd, char *arg) {
  assert(b && arg);
  switch (cmd) {
  case 'N':
    b->n = atol(arg);
    break;
  case 'G':
    b->type = (block_type_t)atoi(arg);
    break;
  case 'X':
    point_set_x(b->target, atof(arg));
    break;
  case 'Y':
    point_set_y(b->target, atof(arg));
    break;
  case 'Z':
    point_set_z(b->target, atof(arg));
    break;
  case 'I':
    b->i = atof(arg);
    break;
  case 'J':
    b->j = atof(arg);
    break;
  case 'R':
    b->r = atof(arg);
    break;
  case 'F': // also support "FMAX"
    if (strcmp(arg, "MAX") == 0) {
      b->feedrate = machine_fmax(b->machine);
    } else {
      b->feedrate = MIN(atof(arg), machine_fmax(b->machine));
    }
    break;
  case 'S':
    b->spindle = atof(arg);
    break;
  case 'T':
    b->tool = atol(arg);
    break;
  default:
    wprintf("Unsupported G-code command \"%c\"\n", cmd);
    return NOCOMMAND_ERR;
  }
  if (b->r && (b->i || b->j)) {
    wprintf("Cannot mix R with I and J\n");
    return ARC_ERR;
  }
  return NO_ERR;
}

static data_t quantize(data_t t, data_t tq, data_t *dq) {
  data_t q;
  q = ((size_t)(t / tq) + 1) * tq;
  *dq = q - t;
  return q;
}

/***************************************************************************************************/
// LOOK_AHEAD Implementation

void block_velocities(block_t *b) {
  assert(b);
  data_t v_i, v_m, v_f;
  data_t alpha, alpha_i;

  data_t seg_a, seg_b, seg_c;

  point_t *vels = point_new();
  if (!vels) {
    eprintf("Could not allocate memory for vels point\n");
  }

  v_i = b->prev ? b->prev->prof->v_f : 0.0;
  v_m = b->arc_feedrate / 60.0;

  
  if (v_m == 0.0) {
    v_f = 0.0;
  }  else if (b->next && block_type(b->next)) {
    // delta of current segment and next segment
    seg_a = point_dist(start_point(b), b->target);
    seg_b = point_dist(start_point(b->next), b->next->target);
    seg_c = point_dist(start_point(b), b->next->target);

    alpha_i = M_PI - acos((pow(seg_a, 2) + pow(seg_b, 2) - pow(seg_c, 2)) / (2 * seg_a * seg_b));
    alpha_i = alpha_i < M_PI_2 ? alpha_i : M_PI_2;
    alpha = MAX(0.0, cos(2 * alpha_i));

    // final velocity - based on scaled average of v_i of current and next segment
    v_f = (v_m + (b->next->arc_feedrate / 60.0)) / 2 * alpha;
  } else {
    v_f = 0.0;
  }

  point_set_xyz(vels, v_i, v_m, v_f);
  char *vel = NULL;
  point_inspect(vels, &vel);
  fprintf(stderr, "%s\n", vel);

  b->prof->v_i = v_i;
  b->prof->v_m = v_m;
  b->prof->v_f = v_f;
}

void block_acceleration(block_t *b) {
  assert(b);

  data_t A;
  data_t ds_1, ds_m, ds_2, s_f; // abscissa
  data_t v_i, v_m, v_f;

  point_t *accl = point_new();

  A = b->acc;
  v_i = b->prof->v_i;
  v_m = b->prof->v_m;
  v_f = b->prof->v_f;
  s_f = b->length;              // total length of segment

  // We have two cases => LONG SEGMENT (v_i < v_m) & SHORT SEGMENT (v)
  //
  // test for short segment
  ds_1 = (pow(v_m, 2) - pow(v_i, 2)) / (2 * A);
  ds_2 = (pow(v_f, 2) - pow(v_m, 2)) / (2 * A);
  if ((ds_1 > s_f) || (ds_2 > s_f)) { // SHORT SEGMENT
    ds_1 = s_f;
    ds_2 = 0;

    // b->prof->v_f = sqrt((2 * A * ds_1) + pow(v_i, 2));

    // if (b->next)
    //   b->next->prof->v_i = b->prof->v_f;

    // iprintf("SHORT SEGMENT\n");
  } else if (v_i < v_m) { // LONG SEGMENT -> [v_i < v_m]
    if (v_f > v_m) {            // A -> v_f > v_m, s_0 < s_1 < s_2 < s_f
      ds_1 = (pow(v_m, 2) - pow(v_i, 2)) / (2 * A);
      ds_2 = (pow(v_f, 2) - pow(v_m, 2)) / (2 * A);
    } else {                    // B -> v_f < v_m, s_0 < s_1 < s_f < s_2
      ds_1 = (pow(v_m, 2) - pow(v_i, 2)) / (2 * A);
      ds_2 = 0;
    }
  // iprintf("LONG SEGMENT\n");
  } else { // LONG SEGMENT -> [v_i > v_m]
    if (v_f > v_m) {            // A -> v_f > v_m, s_1 < s_0 < s_2 < s_f
      ds_1 = 0;
      ds_2 = (pow(v_f, 2) - pow(v_m, 2)) / (2 * A);
    } else {                    // B -> v_f < v_m, s_1 < s_0 < s_f < s_2
      ds_1 = 0;
      ds_2 = 0;
    }

  // iprintf("LONG SEGMENT\n");
  }

  ds_m = s_f - (ds_1 + ds_2);


  point_set_xyz(accl, ds_1, ds_m, ds_2);
  char *acc = NULL;
  point_inspect(accl, &acc);
  // fprintf(stderr, "%s => %f\n", acc, s_f);

  b->prof->ds_1 = ds_1;
  b->prof->ds_2 = ds_2;
  b->prof->ds_m = ds_m;
}

void block_deceleration(block_t *b) {
  assert(b);

  data_t A;
  data_t ds_1, ds_m, ds_2, s_f; // abscissa
  data_t v_i, v_m, v_f;

  point_t *accl = point_new();
  point_t *vels = point_new();
  if (!vels) {
    eprintf("Could not allocate memory for vels point\n");
  }

  // if (!b->prev->type) b


  A = b->acc;
  v_i = b->prof->v_i;
  v_m = b->prof->v_m;           // maintain/mid-level velocity
  v_f = b->prof->v_f;
  s_f = b->length;              // total length of segment

  // ( 9 TOTAL POSSIBLE CASES )
  // POS => 1   |   2   |   3
  //        A       M       A
  //        A                 
  //                M   
  //        D       M       D
  //                        D
  //        A       M       D
  //        A               D
  //        D               A
  //        D       M       A
  // We have two cases => LONG SEGMENT & SHORT SEGMENT
  //
  // test for short segment
  ds_2 = (pow(v_f, 2) - pow(v_m, 2)) / (2 * A);
  ds_m = 0;
  ds_1 = (pow(v_m, 2) - pow(v_i, 2)) / (2 * A);

  if (((fabs(ds_2) > s_f) && (v_f < v_m)) || ((fabs(ds_1) > s_f) && (v_m < v_i))) { // SHORT SEGMENT \\ [D]
    ds_1 = 0;
    ds_2 = -s_f;
    ds_m = 0;

    // last segment/block in program has zero final velocity
    if (!(block_type(b->next)))
      v_f = 0;

    // v_m = 0;
    v_i = sqrt(pow(v_f, 2) - (2 * A * ds_2));
    // if (b->prev)
    //   b->prev->prof->v_f = v_i;

    // iprintf("SHORT SEGMENT -> %f\n", b->prev->prof->v_f);
  } else if (v_f < v_m) { // LONG SEGMENT -> [v_f < v_m] // POSSIBLE CASE ends in [D]
    // v_m  =  sqrt(pow(v_f, 2) - (2 * A * ds_2));

    if (v_i > v_m) {            // A -> v_i > v_m, s_0 < s_1 < s_2 < s_f               [DMD / MD]
      ds_1 = (pow(v_m, 2) - pow(v_i, 2)) / (2 * A);
      ds_m = s_f - (- ds_1 - ds_2);
      if (ds_m < 0) {
        ds_1 = 0;
        ds_m = 0;
        ds_2 = -s_f;
        v_i = sqrt(pow(v_f, 2) - (2 * A * ds_2));
      }
    } else if (v_i < v_m) {                    // B -> v_i < v_m, s_0 < s_1 < s_f < s_2   [AMD / AD]
      ds_1 = b->prof->ds_1;
      ds_m = s_f - (ds_1 - ds_2);
      if (ds_m < 0) {
        ds_1 = s_f / 2;
        ds_m = 0;
        ds_2 = -s_f / 2;
        // v_m = sqrt(pow(v_f, 2) - (2 * A * ds_2));
        v_i = sqrt(pow(v_m, 2) - (2 * A * ds_1));
      }
    } 
    else {
      ds_1 = 0;
      ds_m = s_f - (ds_1 - ds_2);
      // b->prev->prof->v_f = 
    }
  // iprintf("LONG SEGMENT\n");
  } else if (v_f > v_m) { // LONG SEGMENT -> [v_f > v_m] // POSSIBLE CASE ends in [A]
    if (v_i > v_m) {            // A -> v_i > v_m, s_1 < s_0 < s_2 < s_f                  [DMA / DA]
      ds_1 = (pow(v_m, 2) - pow(v_i, 2)) / (2 * A);
      // ds_2 = b->prof->ds_2;
      ds_m = s_f - (- ds_1 + ds_2);
      if (ds_m < 0) {
        ds_1 = -s_f / 2;
        ds_m = 0;
        ds_2 = s_f / 2;
        // v_m = sqrt(pow(v_f, 2) - (2 * A * ds_2));
        v_i = sqrt(pow(v_m, 2) - (2 * A * ds_1));
      }
    } else if (v_i < v_m) {                    // B -> v_i < v_m, s_1 < s_0 < s_f < s_2
      ds_1 = b->prof->ds_1;
      // ds_2 = b->prof->ds_2;
      ds_m = s_f - (ds_1 + ds_2);;
      if (ds_m < 0) {
        ds_m = 0;
        ds_1 = 0;
        ds_2 = s_f;
        v_i = sqrt(pow(v_f, 2) - (2 * A * ds_2));
        // b->prof->v_i = v_i;
        // block_acceleration(b);
        // ds_1 = b->prof->ds_1;
        // ds_m = b->prof->ds_m;
        // ds_2 = b->prof->ds_2;
        // v_m = sqrt(pow(v_f, 2) - (2 * A * ds_2));
        // v_i = sqrt(pow(v_m, 2) - (2 * A * ds_1));
      }
      // else {
      //   ds_1 = b->prof->ds_1;
      //   ds_m = s_f - (ds_1 + ds_2);
      //   if (ds_m < 0) {
      //     ds_m = 0;
      //     ds_1 = 0;
      //     ds_2 = s_f;
      //     v_i = sqrt(pow(v_f, 2) - (2 * A * ds_2));
      //   }
      // }
    }
  // iprintf("LONG SEGMENT\n");
  } else {
    // ds_1 = b->prof->ds_1;
    // ds_2 = b.pf
    // ds_m = s_f - (ds_1 + ds_2);;
    // if (ds_m < 0) {
    //   ds_1 = s_f;
    //   v_i = sqrt(pow(v_f, 2) - (2 * A * ds_1));
    //   b->prof->v_i = v_i;
    //   block_acceleration(b);
      ds_1 = b->prof->ds_1;
      ds_m = b->prof->ds_m;
      ds_2 = b->prof->ds_2;
    // }
  }

  if (b->prev)
    b->prev->prof->v_f = v_i;

  point_set_xyz(accl, ds_1, ds_m, ds_2);
  char *acc = NULL;
  point_inspect(accl, &acc);
  fprintf(stderr, "%s => %f\n", acc, s_f);

  point_set_xyz(vels, v_i, v_m, v_f);
  char *vel = NULL;
  point_inspect(vels, &vel);
  fprintf(stderr, "%s\n\n", vel);


  b->prof->ds_1 = ds_1;
  b->prof->ds_2 = ds_2;
  b->prof->ds_m = ds_m;

  b->prof->v_i = v_i;
  b->prof->v_m = v_m;
  b->prof->v_f = v_f;

}

void block_abs_to_time(block_t *b) {
  assert(b);

  data_t A;
  data_t ds_1, ds_m, ds_2, s_f; // abscissa
  data_t dt, dt_1, dt_2, dt_m, dq;
  data_t v_i, v_m, v_f;


  A = b->acc;
  v_i = b->prof->v_i;
  v_m = b->prof->v_m;           // maintain/mid-level velocity
  v_f = b->prof->v_f;
  s_f = b->length;              // total length of segment
  ds_1 = b->prof->ds_1;
  ds_m = b->prof->ds_m;
  ds_2 = b->prof->ds_2;

  dt_m = 0;

  dt_1 = curvlin_to_time(v_i, A, ds_1);
  printf("\ndt1 of this block is => %.04f | Length => %.03f | v_i => %.02f\n", dt_1, s_f, v_i);
  
  if (ds_2 != 0) {
    if (ds_m > 0) {
      dt_2 = ds_2 > 0 ? curvlin_to_time(v_m, A, ds_2) : curvlin_to_time(v_m, -A, ds_2);
      dt_m = s_f / v_m - (dt_1 + dt_2) / 2.0;
    }
    else
      dt_2 = ds_2 > 0 ? curvlin_to_time(v_i, A, ds_2) : curvlin_to_time(v_i, -A, ds_2);
  }
  else
    dt_2 = 0.0;

  printf("dtm of this block is => %.04f | Length => %.03f | v_m => %.02f\n", dt_m, s_f, v_m);
  printf("dt2 of this block is => %.04f | Length => %.03f | v_f => %.02f\n", dt_2, s_f, v_f);

  b->prof->dt_1 = dt_1;
  b->prof->dt_2 = dt_2;
  b->prof->dt_m = dt_m;
//   b->prof->a = a;
//   b->prof->d = d;
//   b->prof->f = f_m;
//   b->prof->dt = dt;
//   b->prof->l = l;

}

static data_t curvlin_to_time(data_t v_0, data_t a, data_t ds) {
  data_t t1, t2;
  t1 = (1 / a) * (-v_0 + sqrt((2 * a * fabs(ds)) + pow(v_0, 2)));
  t2 = (1 / a) * (-v_0 - sqrt((2 * a * fabs(ds)) + pow(v_0, 2)));

  if (t1 > 0) return t1;
  else if (t2 > 0) return t2;
  else return 0;
}
/***************************************************************************************************/


// static void block_compute(block_t *b) {
//   assert(b);
//   data_t A, a, d;
//   data_t dt, dt_1, dt_2, dt_m, dq;
//   data_t f_m, l;

//   A = b->acc;
//   f_m = b->arc_feedrate / 60.0;
//   l = b->length;
//   dt_1 = f_m / A;
//   dt_2 = dt_1;
//   dt_m = l / f_m - (dt_1 + dt_2) / 2.0;

//   if (dt_m > 0) { // Trapezoidal profile
//     dt = quantize(dt_1 + dt_m + dt_2, machine_tq(b->machine), &dq);
//     dt_m = dt_m + dq;
//     f_m = (2 * l) / (dt_1 + dt_2 + 2 * dt_m);
//   } else { // Triangular profile
//     dt_1 = dt_2 = sqrt(l / A);
//     dt = quantize(dt_1 + dt_2, machine_tq(b->machine), &dq);
//     dt_m = 0;
//     dt_2 = dt_2 + dq;
//     f_m = 2 * l / (dt_1 + dt_2);
//   }
//   a = f_m / dt_1;
//   d = -(f_m / dt_2);
//   b->prof->dt_1 = dt_1;
//   b->prof->dt_2 = dt_2;
//   b->prof->dt_m = dt_m;
//   b->prof->a = a;
//   b->prof->d = d;
//   b->prof->f = f_m;
//   b->prof->dt = dt;
//   b->prof->l = l;
// }

static ccnc_error_t block_arc(block_t *b) {

  data_t x0, y0, z0, xc, yc, xf, yf, zf, r;
  point_t *p0 = start_point(b);
  x0 = point_x(p0);
  y0 = point_y(p0);
  z0 = point_z(p0);
  xf = point_x(b->target);
  yf = point_y(b->target);
  zf = point_z(b->target);

  if (b->r) { // if the radius is given
    data_t dx = point_x(b->delta);
    data_t dy = point_y(b->delta);
    r = b->r;
    data_t dxy2 = pow(dx, 2) + pow(dy, 2);
    data_t sq = sqrt(-pow(dy, 2) * dxy2 * (dxy2 - 4 * r * r));
    // signs table
    // sign(r) | CW(-1) | CCW(+1)
    // --------------------------
    //      -1 |     +  |    -
    //      +1 |     -  |    +
    int s = (r > 0) - (r < 0);
    s *= (b->type == CCWA ? 1 : -1);
    xc = x0 + (dx - s * sq / dxy2) / 2.0;
    yc = y0 + dy / 2.0 + s * (dx * sq) / (2 * dy * dxy2);
  } else { // if I,J are given
    data_t r2;
    r = hypot(b->i, b->j);
    xc = x0 + b->i;
    yc = y0 + b->j;
    r2 = hypot(xf - xc, yf - yc);
    if (fabs(r - r2) > machine_error(b->machine)) {
      fprintf(stderr, "Arc endpoints mismatch error (%f)\n", r - r2);
      block_print(b, stderr);
      return 1;
    }
    b->r = r;
  }
  point_set_x(b->center, xc);
  point_set_y(b->center, yc);
  b->theta_0 = atan2(y0 - yc, x0 - xc);
  b->dtheta = atan2(yf - yc, xf - xc) - b->theta_0;
  // we need the net angle so we take the 2PI complement if negative
  if (b->dtheta < 0)
    b->dtheta = 2 * M_PI + b->dtheta;
  // if CW, take the negative complement
  if (b->type == CWA)
    b->dtheta = -(2 * M_PI - b->dtheta);
  //
  b->length = hypot(zf - z0, b->dtheta * b->r);
  // from now on, it's safer to drop the sign of radius angle
  b->r = fabs(b->r);
  return 0;
}



/*
  ____             _                      _       
 | __ )  ___   ___| | __  _ __ ___   __ _(_)_ __  
 |  _ \ / _ \ / __| |/ / | '_ ` _ \ / _` | | '_ \ 
 | |_) | (_) | (__|   <  | | | | | | (_| | | | | |
 |____/ \___/ \___|_|\_\ |_| |_| |_|\__,_|_|_| |_|
                                                  
*/

#ifdef BLOCK_MAIN

int main(int argc, char const **argv) {
  machine_t *m = machine_new(argv[1]);
  block_t *b1 = NULL, *b2 = NULL, *b3 = NULL, *b4 = NULL;
  if (!m) {
    eprintf("Could not create the machine object\n");
    exit(EXIT_FAILURE);
  }

  b1 = block_new("N10 G01 X90 Y90 Z100 T3 F1000", NULL, m);
  b2 = block_new("N20 G01 y100 S2000", b1, m);
  b3 = block_new("N30 G01 Y200", b2, m);
  b4 = block_new("N40 G00 x0 y0 z0", b3, m);

  block_print(b1, stderr);
  block_print(b2, stderr);
  block_print(b3, stderr);
  block_print(b4, stderr);

  wprintf("Interpolation of block N20 (duration: %f s)\n", block_dt(b2));
  {
    data_t t = 0, tq = machine_tq(m), dt = block_dt(b2);
    data_t lambda = 0, v = 0;
    printf("t lambda v x y z\n");
    for (t = 0; t - dt < tq/10.0; t += tq) {
      block_interpolate_t(b2, t, &lambda, &v);
      printf("%f %f %f %f %f %f\n", t, lambda, v, 
        point_x(machine_setpoint(m)),
        point_y(machine_setpoint(m)),
        point_z(machine_setpoint(m))
      );
    }
  }

  block_free(b1);
  block_free(b2);
  block_free(b3);
  block_free(b4);
  machine_free(m);
  return 0;
}

#endif // BLOCK_MAIN