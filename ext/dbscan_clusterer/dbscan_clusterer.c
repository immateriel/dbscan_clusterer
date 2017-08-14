/* Copyright 2015 Gagarine Yaikhom (MIT License) */
/* Modified by Julien Boulnois 2017 */

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define UNCLASSIFIED -1
#define NOISE -2

#define CORE_POINT 1
#define NOT_CORE_POINT 0

#define SUCCESS 0
#define FAILURE -3

#define MAX_POINT_ELEMENTS 8

typedef struct dist_s dist_t;
struct dist_s {
    void *dist;
};

typedef struct point_s point_t;
struct point_s {
    double elements[MAX_POINT_ELEMENTS];
    unsigned char num_elements;
    int cluster_id;
};

typedef struct node_s node_t;
struct node_s {
    unsigned int index;
    node_t *next;
};

typedef struct epsilon_neighbours_s epsilon_neighbours_t;
struct epsilon_neighbours_s {
    unsigned int num_members;
    node_t *head, *tail;
};

node_t *create_node(unsigned int index);
int append_at_end(
     unsigned int index,
     epsilon_neighbours_t *en);
epsilon_neighbours_t *get_epsilon_neighbours(
    unsigned int index,
    point_t *points,
    unsigned int num_points,
    double epsilon,
    double (*dist)(point_t *a, point_t *b));
void print_epsilon_neighbours(
    point_t *points,
    epsilon_neighbours_t *en);
void destroy_epsilon_neighbours(epsilon_neighbours_t *en);
void dbscan(
    point_t *points,
    unsigned int num_points,
    double epsilon,
    unsigned int minpts,
    double (*dist)(point_t *a, point_t *b));
int expand(
    unsigned int index,
    unsigned int cluster_id,
    point_t *points,
    unsigned int num_points,
    double epsilon,
    unsigned int minpts,
    double (*dist)(point_t *a, point_t *b));
int spread(
    unsigned int index,
    epsilon_neighbours_t *seeds,
    unsigned int cluster_id,
    point_t *points,
    unsigned int num_points,
    double epsilon,
    unsigned int minpts,
    double (*dist)(point_t *a, point_t *b));
double euclidean_dist(point_t *a, point_t *b);
double adjacent_intensity_dist(point_t *a, point_t *b);

node_t *create_node(unsigned int index)
{
    node_t *n = (node_t *) calloc(1, sizeof(node_t));
    if (n == NULL)
        perror("Failed to allocate node.");
    else {
        n->index = index;
        n->next = NULL;
    }
    return n;
}

int append_at_end(
     unsigned int index,
     epsilon_neighbours_t *en)
{
    node_t *n = create_node(index);
    if (n == NULL) {
        free(en);
        return FAILURE;
    }
    if (en->head == NULL) {
        en->head = n;
        en->tail = n;
    } else {
        en->tail->next = n;
        en->tail = n;
    }
    ++(en->num_members);
    return SUCCESS;
}

epsilon_neighbours_t *get_epsilon_neighbours(
    unsigned int index,
    point_t *points,
    unsigned int num_points,
    double epsilon,
    double (*dist)(point_t *a, point_t *b))
{
    unsigned int i;
    epsilon_neighbours_t *en = (epsilon_neighbours_t *)
        calloc(1, sizeof(epsilon_neighbours_t));
    if (en == NULL) {
        perror("Failed to allocate epsilon neighbours.");
        return en;
    }
    for (i = 0; i < num_points; ++i) {
        if (i == index)
            continue;
        if (dist(&points[index], &points[i]) > epsilon)
            continue;
        else {
            if (append_at_end(i, en) == FAILURE) {
                destroy_epsilon_neighbours(en);
                en = NULL;
                break;
            }
        }
    }
    return en;
}

void destroy_epsilon_neighbours(epsilon_neighbours_t *en)
{
    if (en) {
        node_t *t, *h = en->head;
        while (h) {
            t = h->next;
            free(h);
            h = t;
        }
        free(en);
    }
}

void dbscan(
    point_t *points,
    unsigned int num_points,
    double epsilon,
    unsigned int minpts,
    double (*dist)(point_t *a, point_t *b))
{
    unsigned int i, cluster_id = 0;
    for (i = 0; i < num_points; ++i) {
        if (points[i].cluster_id == UNCLASSIFIED) {
            if (expand(i, cluster_id, points,
                       num_points, epsilon, minpts,
                       dist) == CORE_POINT)
                ++cluster_id;
        }
    }
}

int expand(
    unsigned int index,
    unsigned int cluster_id,
    point_t *points,
    unsigned int num_points,
    double epsilon,
    unsigned int minpts,
    double (*dist)(point_t *a, point_t *b))
{
    int return_value = NOT_CORE_POINT;
    epsilon_neighbours_t *seeds =
        get_epsilon_neighbours(index, points,
                               num_points, epsilon,
                               dist);
    if (seeds == NULL)
        return FAILURE;

    if (seeds->num_members < minpts)
        points[index].cluster_id = NOISE;
    else {
        points[index].cluster_id = cluster_id;
        node_t *h = seeds->head;
        while (h) {
            points[h->index].cluster_id = cluster_id;
            h = h->next;
        }

        h = seeds->head;
        while (h) {
            spread(h->index, seeds, cluster_id, points,
                   num_points, epsilon, minpts, dist);
            h = h->next;
        }

        return_value = CORE_POINT;
    }
    destroy_epsilon_neighbours(seeds);
    return return_value;
}

int spread(
    unsigned int index,
    epsilon_neighbours_t *seeds,
    unsigned int cluster_id,
    point_t *points,
    unsigned int num_points,
    double epsilon,
    unsigned int minpts,
    double (*dist)(point_t *a, point_t *b))
{
    epsilon_neighbours_t *spread =
        get_epsilon_neighbours(index, points,
                       num_points, epsilon,
                       dist);
    if (spread == NULL)
        return FAILURE;
    if (spread->num_members >= minpts) {
        node_t *n = spread->head;
        point_t *d;
        while (n) {
            d = &points[n->index];
            if (d->cluster_id == NOISE ||
                d->cluster_id == UNCLASSIFIED) {
                if (d->cluster_id == UNCLASSIFIED) {
                    if (append_at_end(n->index, seeds)
                        == FAILURE) {
                        destroy_epsilon_neighbours(spread);
                        return FAILURE;
                    }
                }
                d->cluster_id = cluster_id;
            }
            n = n->next;
        }
    }

    destroy_epsilon_neighbours(spread);
    return SUCCESS;
}

#define d_max(x,y) ((x) >= (y)) ? (x) : (y)

double approximated2d_dist(point_t *a, point_t *b) {
    double dx=abs(a->elements[0] - b->elements[0]);
    double dy=abs(a->elements[1] - b->elements[1]);
    return 0.394 * (dx + dy) + 0.554 * d_max(dx, dy);
}

double euclidean2d_dist(point_t *a, point_t *b) {
    return sqrt(pow(a->elements[0] - b->elements[0], 2) +
    pow(a->elements[1] - b->elements[1], 2));  
}

double euclidean_dist(point_t *a, point_t *b)
{
    int i=0;
    float sum=0.0;
    for(i=0;i<a->num_elements;i++) {
        sum+=pow(a->elements[i] - b->elements[i],2);
    }

    return sqrt(sum);
}

/* ruby binding */

#include <ruby.h>

VALUE DbscanClusterer = Qnil;
void Init_dbscan_clusterer();

//VALUE method_dbscan_clusterer_dbscan(VALUE self, VALUE points, VALUE epsilon, VALUE minpts);
VALUE method_dbscan_clusterer_dbscan(int argc, VALUE *argv, VALUE self);

void Init_dbscan_clusterer() {
    DbscanClusterer = rb_define_module("DbscanClusterer");
    rb_define_singleton_method(DbscanClusterer, "dbscan", method_dbscan_clusterer_dbscan, -1);
}

static inline void _ruby_array_to_point(VALUE rel, point_t *p) {
    int j;
    unsigned char num_elts = RARRAY_LEN(rel);
    for (j = 0; j < num_elts; j++)
        p->elements[j] = NUM2DBL(rb_ary_entry(rel, j));

    p->num_elements = num_elts;
    p->cluster_id = UNCLASSIFIED;

    return p;
}

static inline VALUE _point_to_ruby_array(point_t* p) {
    int j;
    VALUE rpoint = rb_ary_new2(p->num_elements);
    for (j = 0; j < p->num_elements; j++)
        rb_ary_store(rpoint, j, DBL2NUM(p->elements[j]));

    return rpoint;
}
/*
int ruby_dist(point_t *a, point_t *b) {

     VALUE passthrough = (VALUE)userdata;
     VALUE cb;
     VALUE cbdata;

     cb = rb_ary_entry(passthrough, 0);
     cbdata = rb_ary_entry(passthrough, 1);

     rb_funcall(cb, rb_intern("call"), 2, INT2NUM(event),
                cbdata);

     return 0;
}

VALUE ext_register(VALUE obj, VALUE cb) {
     VALUE passthrough;

     if (rb_class_of(cb) != rb_cProc)
         rb_raise(rb_eTypeError, "Expected Proc callback");

     passthrough = rb_ary_new();
     rb_ary_store(passthrough, 0, cb);

     register_async_callback(internal_callback,
                             (void *)passthrough);
 }
*/
// rb_define_method(c_extension, "register_async_proc",
//                  ext_register, 2);
// WARN: there is no check of points elements size
VALUE method_dbscan_clusterer_dbscan(int argc, VALUE *argv, VALUE self) {
//VALUE method_dbscan_clusterer_dbscan(VALUE self, VALUE points, VALUE epsilon, VALUE minpts) {
    unsigned int i,j;
    unsigned int num_points;
    VALUE results;
    point_t *p;
    void *c_dist_func;

    VALUE points;
    VALUE epsilon;
    VALUE minpts;
    VALUE dist_method;
    VALUE dist_proc;

    VALUE euclidean_sym = ID2SYM(rb_intern("euclidean"));
    VALUE euclidean2d_sym = ID2SYM(rb_intern("euclidean2d"));
    VALUE approximated2d_sym = ID2SYM(rb_intern("approximated2d"));
//    VALUE ruby_sym = ID2SYM(rb_intern("ruby"));

    rb_scan_args(argc, argv, "32", &points, &epsilon, &minpts, &dist_method, &dist_proc);

    if(dist_method == Qnil) {
        c_dist_func = euclidean_dist;
    } else if(dist_method == euclidean_sym) {
        c_dist_func = euclidean_dist;        
    } else if(dist_method == euclidean2d_sym) {
        c_dist_func = euclidean2d_dist;
    } else if(dist_method == approximated2d_sym) {
        c_dist_func = approximated2d_dist;
    } else {
        c_dist_func = euclidean_dist;        
    }

    num_points = RARRAY_LEN(points);
    p = (point_t *) calloc(num_points, sizeof(point_t));

    for(i=0;i<num_points;i++) {
        VALUE rel = rb_ary_entry(points, i);
        _ruby_array_to_point(rel,&p[i]);
    }

    dbscan(p, num_points, NUM2DBL(epsilon), NUM2INT(minpts), c_dist_func);

    results = rb_hash_new();

    for(i=0;i<num_points;i++) {
        VALUE rpoint = _point_to_ruby_array(&p[i]);
        VALUE rid = INT2NUM(p[i].cluster_id);

        VALUE rpoints = rb_hash_aref(results, rid);
        if(NIL_P(rpoints)) {
            rpoints = rb_ary_new();
            rb_hash_aset(results, rid, rpoints);
        }
        rb_ary_push(rpoints, rpoint);
    }

    free(p);

    return results;
}
