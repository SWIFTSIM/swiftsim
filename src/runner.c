/*******************************************************************************
 * This file is part of GadgetSMP.
 * Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 ******************************************************************************/

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <omp.h>
#include <sched.h>

/* Local headers. */
#include "cycle.h"
#include "lock.h"
#include "task.h"
#include "part.h"
#include "cell.h"
#include "space.h"
#include "queue.h"
#include "runner.h"
#include "runner_iact.h"

/* Error macro. */
#define error(s) { printf( "%s:%s:%i: %s\n" , __FILE__ , __FUNCTION__ , __LINE__ , s ); abort(); }

/* Convert cell location to ID. */
#define cell_getid( cdim , i , j , k ) ( (int)(k) + (cdim)[2]*( (int)(j) + (cdim)[1]*(int)(i) ) )

/* The timers. */
ticks runner_timer[ runner_timer_count ];

/* The counters. */
int runner_counter[ runner_counter_count ];

        

const float runner_shift[13*3] = {
     5.773502691896258e-01 ,  5.773502691896258e-01 ,  5.773502691896258e-01 ,
     7.071067811865475e-01 ,  7.071067811865475e-01 ,  0.0                   ,
     5.773502691896258e-01 ,  5.773502691896258e-01 , -5.773502691896258e-01 ,
     7.071067811865475e-01 ,  0.0                   ,  7.071067811865475e-01 ,
     1.0                   ,  0.0                   ,  0.0                   ,
     7.071067811865475e-01 ,  0.0                   , -7.071067811865475e-01 ,
     5.773502691896258e-01 , -5.773502691896258e-01 ,  5.773502691896258e-01 ,
     7.071067811865475e-01 , -7.071067811865475e-01 ,  0.0                   ,
     5.773502691896258e-01 , -5.773502691896258e-01 , -5.773502691896258e-01 ,
     0.0                   ,  7.071067811865475e-01 ,  7.071067811865475e-01 ,
     0.0                   ,  1.0                   ,  0.0                   ,
     0.0                   ,  7.071067811865475e-01 , -7.071067811865475e-01 ,
     0.0                   ,  0.0                   ,  1.0                   ,
    };
const char runner_flip[27] = { 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 0 ,
                               0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }; 


/* Import the density functions. */
#define FUNCTION density
#include "runner_doiact.h"

#undef FUNCTION
#define FUNCTION force
#include "runner_doiact.h"


/** 
 * @brief Sort the tasks in topological order over all queues.
 *
 * @param r The #runner.
 */
 
void runner_ranktasks ( struct runner *r ) {

    int i, j = 0, k, temp, left = 0, rank;
    struct task *t;
    struct space *s = r->s;
    int *tid;

    /* Run throught the tasks and get all the waits right. */
    for ( k = 0 ; k < s->nr_tasks ; k++ ) {
        for ( j = 0 ; j < s->tasks[k].nr_unlock_tasks ; j++ )
            s->tasks[k].unlock_tasks[j]->wait += 1;
        }
        
    /* Allocate and init the task-ID array. */
    if ( ( tid = (int *)malloc( sizeof(int) * s->nr_tasks ) ) == NULL )
        error( "Failed to allocate temporary tid array." );
    for ( k = 0 ; k < s->nr_tasks ; k++ )
        tid[k] = k;
        
    /* Main loop. */
    for ( rank = 0 ; left < s->nr_tasks ; rank++ ) {
        
        /* Load the tids of tasks with no waits. */
        for ( k = left ; k < s->nr_tasks ; k++ )
            if ( s->tasks[ tid[k] ].wait == 0 ) {
                temp = tid[j]; tid[j] = tid[k]; tid[k] = temp;
                j += 1;
                }

        /* Traverse the task tree and add tasks with no weight. */
        for ( i = left ; i < j ; i++ ) {
            t = &s->tasks[ tid[i] ];
            t->rank = rank;
            s->tasks_ind[i] = t - s->tasks;
            /* printf( "runner_ranktasks: task %i of type %s has rank %i.\n" , i , 
                (t->type == task_type_self) ? "self" : (t->type == task_type_pair) ? "pair" : "sort" , rank ); */
            for ( k = 0 ; k < t->nr_unlock_tasks ; k++ )
                t->unlock_tasks[k]->wait -= 1;
            }
            
        /* The new left (no, not tony). */
        left = j;
            
        }
        
    /* Release the temporary array. */
    free(tid);
    
    }


/**
 * @brief Sort the entries in ascending order using QuickSort.
 *
 * @param sort The entries
 * @param N The number of entries.
 */
 
void runner_dosort_ascending ( struct entry *sort , int N ) {

    struct {
        short int lo, hi;
        } qstack[10];
    int qpos, i, j, lo, hi, imin;
    struct entry temp;
    float pivot;
        
    /* Sort parts in cell_i in decreasing order with quicksort */
    qstack[0].lo = 0; qstack[0].hi = N - 1; qpos = 0;
    while ( qpos >= 0 ) {
        lo = qstack[qpos].lo; hi = qstack[qpos].hi;
        qpos -= 1;
        if ( hi - lo < 15 ) {
            for ( i = lo ; i < hi ; i++ ) {
                imin = i;
                for ( j = i+1 ; j <= hi ; j++ )
                    if ( sort[j].d < sort[imin].d )
                        imin = j;
                if ( imin != i ) {
                    temp = sort[imin]; sort[imin] = sort[i]; sort[i] = temp;
                    }
                }
            }
        else {
            pivot = sort[ ( lo + hi ) / 2 ].d;
            i = lo; j = hi;
            while ( i <= j ) {
                while ( sort[i].d < pivot ) i++;
                while ( sort[j].d > pivot ) j--;
                if ( i <= j ) {
                    if ( i < j ) {
                        temp = sort[i]; sort[i] = sort[j]; sort[j] = temp;
                        }
                    i += 1; j -= 1;
                    }
                }
            if ( j > ( lo + hi ) / 2 ) {
                if ( lo < j ) {
                    qpos += 1;
                    qstack[qpos].lo = lo;
                    qstack[qpos].hi = j;
                    }
                if ( i < hi ) {
                    qpos += 1;
                    qstack[qpos].lo = i;
                    qstack[qpos].hi = hi;
                    }
                }
            else {
                if ( i < hi ) {
                    qpos += 1;
                    qstack[qpos].lo = i;
                    qstack[qpos].hi = hi;
                    }
                if ( lo < j ) {
                    qpos += 1;
                    qstack[qpos].lo = lo;
                    qstack[qpos].hi = j;
                    }
                }
            }
        }
                
    }
    
    
/**
 * @brief inline helper fuction to merge two entry arrays (forward).
 *
 * @param one the first array
 * @param none the length of the first array
 * @param two the second array
 * @param ntwo the length of the second array
 * @param dest the destination array.
 */
 
inline void merge_forward ( struct entry *__restrict__ one , int none , struct entry *__restrict__ two , int ntwo , struct entry *__restrict__ dest ) {

    int i = 0, j = 0, k = 0;
    
    while ( j < none && k < ntwo )
        if ( one[j].d < two[k].d )
            dest[i++] = one[j++];
        else
            dest[i++] = two[k++];
    if ( j == none )
        for ( ; k < ntwo ; k++ )
            dest[i++] = two[k];
    else
        for ( ; j < none ; j++ )
            dest[i++] = one[j];

    }
    

/**
 * @brief inline helper fuction to merge two entry arrays (forward).
 *
 * @param one the first array
 * @param none the length of the first array
 * @param two the second array
 * @param ntwo the length of the second array
 * @param dest the destination array.
 */
 
inline void merge_backward ( struct entry *__restrict__ one , int none , struct entry *__restrict__ two , int ntwo , struct entry *__restrict__ dest ) {

    int i = none + ntwo - 1, j = none - 1, k = ntwo - 1;
    
    while ( j >= 0 && k >= 0 )
        if ( one[j].d > two[k].d )
            dest[i--] = one[j--];
        else
            dest[i--] = two[k--];
    if ( j < 0 )
        for ( ; k >= 0 ; k-- )
            dest[i--] = two[k];
    else
        for ( ; j >= 0 ; j-- )
            dest[i--] = one[j];

    }
    

/**
 * @brief Sort the particles in the given cell along all cardinal directions.
 *
 * @param r The #runner.
 * @param c The #cell.
 */
 
void runner_dosort ( struct runner_thread *rt , struct cell *c , int flags ) {

    struct entry *finger;
    struct entry *fingers[8];
    struct part *parts = c->parts;
    int j, k, count = c->count;
    int cone, ctwo;
    int i, ind, off[8], inds[8], temp_i;
    // float shift[3];
    float buff[8], px[3];
    struct cell *temp_c;
    TIMER_TIC
    
    /* Does this cell even need to be sorted? */
    for ( temp_c = c ; temp_c != NULL && temp_c->nr_pairs == 0 ; temp_c = temp_c->parent );
    if ( temp_c == NULL )
        return;

    /* start by allocating the entry arrays. */
    if ( lock_lock( &c->lock ) != 0 )
        error( "Failed to lock cell." );
    if ( c->sort == NULL )
        if ( ( c->sort = (struct entry *)malloc( sizeof(struct entry) * (c->count + 1) * 13 ) ) == NULL )
            error( "Failed to allocate sort memory." );
    if ( lock_unlock( &c->lock ) != 0 )
        error( "Failed to unlock cell." );
        
    /* Does this cell have any progeny? */
    if ( c->split ) {
    
        /* Loop over the 13 different sort arrays. */
        for ( j = 0 ; j < 13 ; j++ ) {
        
            /* Has this sort array been flagged? */
            if ( !( flags & (1 << j) ) )
                continue;
                
            if ( 0 ) {
            
                /* Get a finger on the sorting array. */
                finger = &c->sort[ j*(count + 1) ];
                
                /* Merge the two first sub-cells forward into this cell. */
                cone = c->progeny[0]->count;
                ctwo = c->progeny[1]->count;
                merge_forward( &c->progeny[0]->sort[ j*(cone + 1) ] , cone ,
                               &c->progeny[1]->sort[ j*(ctwo + 1) ] , ctwo ,
                               finger );
                               
                /* Merge-in the remaining arrays, alternating forward and
                   backward merges. */
                for ( k = 2 ; k < 8 ; k++ ) {
                    cone = cone + ctwo;
                    ctwo = c->progeny[k]->count;
                    if ( k & 1 )
                        merge_forward( &finger[ count - cone ] , cone ,
                                       &c->progeny[k]->sort[ j*(ctwo + 1) ] , ctwo ,
                                       finger );
                    else
                        merge_backward( finger , cone ,
                                        &c->progeny[k]->sort[ j*(ctwo + 1) ] , ctwo ,
                                        &finger[ count - cone - ctwo ] );
                    }
                
                }
                
            else {
            
                /* Init the particle index offsets. */
                for ( off[0] = 0 , k = 1 ; k < 8 ; k++ )
                    if ( c->progeny[k-1] != NULL )
                        off[k] = off[k-1] + c->progeny[k-1]->count;
                    else
                        off[k] = off[k-1];

                /* Init the entries and indices. */
                for ( k = 0 ; k < 8 ; k++ ) {
                    inds[k] = k;
                    if ( c->progeny[k] != NULL && c->progeny[k]->count > 0 ) {
                        fingers[k] = &c->progeny[k]->sort[ j*(c->progeny[k]->count + 1) ];
                        buff[k] = fingers[k]->d;
                        off[k] = off[k];
                        }
                    else
                        buff[k] = FLT_MAX;
                    }

                /* Sort the buffer. */
                for ( i = 0 ; i < 7 ; i++ )
                    for ( k = i+1 ; k < 8 ; k++ )
                        if ( buff[ inds[k] ] < buff[ inds[i] ] ) {
                            temp_i = inds[i]; inds[i] = inds[k]; inds[k] = temp_i;
                            }

                /* For each entry in the new sort list. */
                finger = &c->sort[ j*(count + 1) ];
                for ( ind = 0 ; ind < count ; ind++ ) {

                    /* Copy the minimum into the new sort array. */
                    finger[ind].d = buff[inds[0]];
                    finger[ind].i = fingers[inds[0]]->i + off[inds[0]];

                    /* Update the buffer. */
                    fingers[inds[0]] += 1;
                    buff[inds[0]] = fingers[inds[0]]->d;

                    /* Find the smallest entry. */
                    for ( k = 1 ; k < 8 && buff[inds[k]] < buff[inds[k-1]] ; k++ ) {
                        temp_i = inds[k-1]; inds[k-1] = inds[k]; inds[k] = temp_i;
                        }

                    } /* Merge. */
                    
                }
            
            /* Add a sentinel. */
            c->sort[ j*(c->count + 1) + c->count ].d = FLT_MAX;
            c->sort[ j*(c->count + 1) + c->count ].i = 0;
            
            } /* loop over sort arrays. */
    
        } /* progeny? */
        
    /* Otherwise, just sort. */
    // else {
    // 
    //     /* Loop over the different cell axes. */
    //     for ( j = 0 ; j < 13 ; j++ ) {
    //     
    //         /* Has this sort array been flagged? */
    //         if ( !( flags & (1 << j) ) )
    //             continue;
    //     
    //         /* Get the shift vector. */
    //         shift[0] = runner_shift[ 3*j + 0 ];
    //         shift[1] = runner_shift[ 3*j + 1 ];
    //         shift[2] = runner_shift[ 3*j + 2 ];
    //         
    //         /* Fill the sort array. */
    //         finger = &c->sort[ j*(count + 1) ];
    //         for ( k = 0 ; k < count ; k++ ) {
    //             finger[k].i = k;
    //             finger[k].d = parts[k].x[0]*shift[0] + parts[k].x[1]*shift[1] + parts[k].x[2]*shift[2];
    //             }
    //             
    //         /* Add the sentinel. */
    //         finger[ c->count ].d = FLT_MAX;
    //         finger[ c->count ].i = 0;
    //             
    //         /* Sort descending. */
    //         runner_dosort_ascending( finger , c->count );
    //     
    //         }
    //         
    //     }
        
    /* Otherwise, just sort. */
    else {
    
        /* Fill the sort array. */
        for ( k = 0 ; k < count ; k++ ) {
            px[0] = parts[k].x[0];
            px[1] = parts[k].x[1];
            px[2] = parts[k].x[2];
            for ( j = 0 ; j < 13 ; j++ )
                if ( flags & (1 << j) ) {
                    c->sort[ j*(count + 1) + k].i = k;
                    c->sort[ j*(count + 1) + k].d = px[0]*runner_shift[ 3*j + 0 ] + px[1]*runner_shift[ 3*j + 1 ] + px[2]*runner_shift[ 3*j + 2 ];
                    }
            }

        /* Add the sentinel and sort. */
        for ( j = 0 ; j < 13 ; j++ )
            if ( flags & (1 << j) ) {
                c->sort[ j*(count + 1) + c->count ].d = FLT_MAX;
                c->sort[ j*(count + 1) + c->count ].i = 0;
                runner_dosort_ascending( &c->sort[ j*(count + 1) ] , c->count );
                }
            
        }
        
    /* Verify the sorting. */
    /* for ( j = 0 ; j < 13 ; j++ ) {
        if ( !( flags & (1 << j) ) )
            continue;
        finger = &c->sort[ j*(c->count + 1) ];
        for ( k = 1 ; k < c->count ; k++ ) {
            if ( finger[k].d < finger[k-1].d )
                error( "Sorting failed, ascending array." );
            if ( finger[k].i >= c->count )
                error( "Sorting failed, indices borked." );
            }
        } */

    #ifdef TIMER_VERBOSE
        printf( "runner_dosort[%02i]: %i parts at depth %i (flags = %i%i%i%i%i%i%i%i%i%i%i%i%i) took %.3f ms.\n" ,
            rt->id , c->count , c->depth ,
            (flags & 0x1000) >> 12 , (flags & 0x800) >> 11 , (flags & 0x400) >> 10 , (flags & 0x200) >> 9 , (flags & 0x100) >> 8 , (flags & 0x80) >> 7 , (flags & 0x40) >> 6 , (flags & 0x20) >> 5 , (flags & 0x10) >> 4 , (flags & 0x8) >> 3 , (flags & 0x4) >> 2 , (flags & 0x2) >> 1 , (flags & 0x1) >> 0 , 
            ((double)TIMER_TOC(runner_timer_dosort)) / CPU_TPS * 1000 ); fflush(stdout);
    #else
        TIMER_TOC(runner_timer_dosort);
    #endif

    }


/**
 * @brief Implements a barrier for the #runner threads.
 *
 */
 
void runner_barrier( struct runner *r ) {

    /* First, get the barrier mutex. */
    if ( pthread_mutex_lock( &r->barrier_mutex ) != 0 )
        error( "Failed to get barrier mutex." );
        
    /* Wait for the barrier to close. */
    while ( r->barrier_count < 0 )
        if ( pthread_cond_wait( &r->barrier_cond , &r->barrier_mutex ) != 0 )
            error( "Eror waiting for barrier to close." );
        
    /* Once I'm in, increase the barrier count. */
    r->barrier_count += 1;
    
    /* If all threads are in, send a signal... */
    if ( r->barrier_count == r->nr_threads )
        if ( pthread_cond_broadcast( &r->barrier_cond ) != 0 )
            error( "Failed to broadcast barrier full condition." );
        
    /* Wait for barrier to be released. */
    while ( r->barrier_count > 0 )
        if ( pthread_cond_wait( &r->barrier_cond , &r->barrier_mutex ) != 0 )
            error( "Error waiting for barrier to be released." );
            
    /* Decrease the counter before leaving... */
    r->barrier_count += 1;
    
    /* If I'm the last one out, signal the condition again. */
    if ( r->barrier_count == 0 )
        if ( pthread_cond_broadcast( &r->barrier_cond ) != 0 )
            error( "Failed to broadcast empty barrier condition." );
            
    /* Last but not least, release the mutex. */
    if ( pthread_mutex_unlock( &r->barrier_mutex ) != 0 )
        error( "Failed to get unlock the barrier mutex." );

    }
    
    
/**
 * @brief The #runner main thread routine.
 *
 * @param data A pointer to this thread's data.
 */
 
void *runner_main ( void *data ) {

    struct runner_thread *rt = (struct runner_thread *)data;
    struct runner *r = rt->r;
    int threadID = rt->id;
    int k, qid, naq, keep, tpq;
    struct queue *queues[ r->nr_queues ], *myq;
    struct task *t;
    struct cell *ci, *cj;
    unsigned int myseed = rand() + rt->id;
    #ifdef TIMER
        ticks stalled;
    #endif
    
    /* Main loop. */
    while ( 1 ) {
    
        /* Wait at the barrier. */
        runner_barrier( r );
        
        /* Set some convenient local data. */
        keep = r->policy & runner_policy_keep;
        myq = &r->queues[ threadID % r->nr_queues ];
        tpq = ceil( ((double)r->nr_threads) / r->nr_queues );
        stalled = 0;
        
        /* Set up the local list of active queues. */
        naq = r->nr_queues;
        for ( k = 0 ; k < naq ; k++ )
            queues[k] = &r->queues[k];
    
        /* Set up the local list of active queues. */
        naq = r->nr_queues;
        for ( k = 0 ; k < naq ; k++ )
            queues[k] = &r->queues[k];
    
        /* Loop while there are tasks... */
        while ( 1 ) {
        
            /* Remove any inactive queues. */
            for ( k = 0 ; k < naq ; k++ )
                if ( queues[k]->next == queues[k]->count ) {
                    naq -= 1;
                    queues[k] = queues[naq];
                    k -= 1;
                    }
            if ( naq == 0 )
                break;
        
            /* Get a task, how and from where depends on the policy. */
            TIMER_TIC
            t = NULL;
            if ( r->nr_queues == 1 ) {
                t = queue_gettask( &r->queues[0] , 1 , 0 );
                }
            else if ( r->policy & runner_policy_steal ) {
                if ( ( myq->next == myq->count ) ||
                     ( t = queue_gettask_new( myq , rt->id , 0 , 0 ) ) == NULL ) {
                    TIMER_TIC2
                    qid = rand_r( &myseed ) % naq;
                    keep = ( r->policy & runner_policy_keep ) &&
                           ( myq->count <= myq->size-tpq );
                    if ( myq->next == myq->count )
                        COUNT(runner_counter_steal_empty);
                    else
                        COUNT(runner_counter_steal_stall);
                    t = queue_gettask_new( queues[qid] , rt->id , 0 , keep );
                    if ( t != NULL && keep )
                        queue_insert( myq , t );
                    TIMER_TOC2(runner_timer_steal);
                    }
                }
            else if ( r->policy & runner_policy_rand ) {
                qid = rand_r( &myseed ) % naq;
                t = queue_gettask( queues[qid] , r->policy & runner_policy_block , 0 );
                }
            else {
                t = queue_gettask( &r->queues[threadID] , r->policy & runner_policy_block , 0 );
                }
            TIMER_TOC(runner_timer_getpair);
            
            /* Did I get anything? */
            if ( t == NULL ) {
                COUNT(runner_counter_stall);
                if ( !stalled )
                    stalled = getticks();
                continue;
                }
            #ifdef TIMER
            else if ( stalled ) {
                stalled = getticks() - stalled;
                __sync_add_and_fetch( &runner_timer[runner_timer_stalled] , stalled );
                #ifdef TIMER_VERBOSE
                    printf( "runner_main[%02i]: stalled %.3f ms\n" , rt->id , ((double)stalled) / CPU_TPS * 1000 );
                    fflush(stdout);
                #endif
                stalled = 0;
                }
            #endif
        
            /* Get the cells. */
            ci = t->ci;
            cj = t->cj;
            
            /* Different types of tasks... */
            switch ( t->type ) {
                case task_type_self:
                    if ( t->subtype == task_subtype_density )
                        runner_doself_density( rt , ci );
                    else if ( t->subtype == task_subtype_force )
                        runner_doself_force( rt , ci );
                    else
                        error( "Unknown task subtype." );
                    cell_unlocktree( ci );
                    break;
                case task_type_pair:
                    if ( t->subtype == task_subtype_density )
                        runner_dopair_density( rt , ci , cj );
                    else if ( t->subtype == task_subtype_force )
                        runner_dopair_force( rt , ci , cj );
                    else
                        error( "Unknown task subtype." );
                    cell_unlocktree( ci );
                    cell_unlocktree( cj );
                    break;
                case task_type_sort:
                    runner_dosort( rt , ci , t->flags );
                    break;
                case task_type_sub:
                    if ( t->subtype == task_subtype_density )
                        runner_dosub_density( rt , ci , cj , t->flags );
                    else if ( t->subtype == task_subtype_force )
                        runner_dosub_force( rt , ci , cj , t->flags );
                    else
                        error( "Unknown task subtype." );
                    cell_unlocktree( ci );
                    if ( cj != NULL )
                        cell_unlocktree( cj );
                    break;
                case task_type_ghost:
                    break;
                default:
                    error( "Unknown task type." );
                }
                
            t->done = 1;
            
            /* Resolve any dependencies. */
            for ( k = 0 ; k < t->nr_unlock_tasks ; k++ )
                if ( __sync_fetch_and_sub( &t->unlock_tasks[k]->wait , 1 ) == 0 )
                    abort();
            for ( k = 0 ; k < t->nr_unlock_cells ; k++ )
                __sync_fetch_and_sub( &t->unlock_cells[k]->wait , 1 );
        
            } /* main loop. */
            
    	/* Any leftover stalls? */    
        #ifdef TIMER
        if ( stalled ) {
            stalled = getticks() - stalled;
            __sync_add_and_fetch( &runner_timer[runner_timer_stalled] , stalled );
            #ifdef TIMER_VERBOSE
                printf( "runner_main[%02i]: stalled %.3f ms\n" , rt->id , ((double)stalled) / CPU_TPS * 1000 );
                fflush(stdout);
            #endif
            stalled = 0;
            }
        #endif
            
        }
        
    /* Be kind, rewind. */
    return NULL;

    }
    

/**
 * @brief Let the #runner loose on the given #space.
 *
 * @param r The #runner.
 * @param s The #space.
 */
 
void runner_run ( struct runner *r , int sort_queues ) {

    int j, k;
    struct space *s = r->s;
    
    /* Run throught the tasks and get all the waits right. */
    for ( k = 0 ; k < s->nr_tasks ; k++ ) {
        s->tasks[k].done = 0;
        for ( j = 0 ; j < s->tasks[k].nr_unlock_tasks ; j++ )
            s->tasks[k].unlock_tasks[j]->wait += 1;
        for ( j = 0 ; j < s->tasks[k].nr_unlock_cells ; j++ )
            s->tasks[k].unlock_cells[j]->wait += 1;
        }
    
    /* Re-set the queues.*/
    if ( sort_queues ) {
        #pragma omp parallel for default(none), shared(r)
        for ( k = 0 ; k < r->nr_queues ; k++ ) {
            queue_sort( &r->queues[k] );
            r->queues[k].next = 0;
            }
        }
    else
        for ( k = 0 ; k < r->nr_queues ; k++ )
            r->queues[k].next = 0;
    
    /* Cry havoc and let loose the dogs of war. */
    r->barrier_count = -r->barrier_count;
    if ( pthread_cond_broadcast( &r->barrier_cond ) != 0 )
        error( "Failed to broadcast barrier open condition." );
        
    /* Sit back and wait for the runner_threads to come home. */
    while ( r->barrier_count < r->nr_threads )
        if ( pthread_cond_wait( &r->barrier_cond , &r->barrier_mutex ) != 0 )
            error( "Error while waiting for barrier." );
    
    }
    
    
/**
 * @brief init a runner with the given number of threads, queues, and
 *      the given policy.
 *
 * @param r The #runner.
 * @param s The #space in which this #runner will run.
 * @param nr_threads The number of threads to spawn.
 * @param nr_queues The number of task queues to create.
 * @param policy The queueing policy to use.
 */
 
void runner_init ( struct runner *r , struct space *s , int nr_threads , int nr_queues , int policy ) {

    #if defined(HAVE_SETAFFINITY)
        cpu_set_t cpuset;
    #endif
    int k, qid, nrq;
    
    /* Store the values. */
    r->s = s;
    r->nr_threads = nr_threads;
    r->nr_queues = nr_queues;
    r->policy = policy;
    
    /* First of all, init the barrier and lock it. */
    if ( pthread_mutex_init( &r->barrier_mutex , NULL ) != 0 )
        error( "Failed to initialize barrier mutex." );
    if ( pthread_cond_init( &r->barrier_cond , NULL ) != 0 )
        error( "Failed to initialize barrier condition variable." );
    if ( pthread_mutex_lock( &r->barrier_mutex ) != 0 )
        error( "Failed to lock barrier mutex." );
    r->barrier_count = 0;
    
    /* Allocate the queues. */
    if ( posix_memalign( (void *)(&r->queues) , 64 , nr_queues * sizeof(struct queue) ) != 0 )
        error( "Failed to allocate queues." );
    bzero( r->queues , nr_queues * sizeof(struct queue) );
        
    /* Init the queues. */
    for ( k = 0 ; k < nr_queues ; k++ )
        queue_init( &r->queues[k] , s->nr_tasks , s->tasks );
        
    /* Rank the tasks in topological order. */
    runner_ranktasks( r );
    
    /* How many queues to fill initially? */
    for ( nrq = 0 , k = nr_queues ; k > 0 ; k = k / 2 )
        nrq += 1;
        
    /* Fill the queues (round-robin). */
    for ( k = 0 ; k < s->nr_tasks ; k++ ) {
        if ( s->tasks[ s->tasks_ind[k] ].type == task_type_none )
            continue;
        // qid = 0;
        // qid = k % nrq;
        qid = k % nr_queues;
        r->queues[qid].tid[ r->queues[qid].count ] = s->tasks_ind[k];
        r->queues[qid].count += 1;
        }
        
    /* Sort the queues topologically. */
    for ( k = 0 ; k < nr_queues ; k++ )
        queue_sort( &r->queues[k] );
        
    /* Allocate and init the threads. */
    if ( ( r->threads = (struct runner_thread *)malloc( sizeof(struct runner_thread) * nr_threads ) ) == NULL )
        error( "Failed to allocate threads array." );
    for ( k = 0 ; k < nr_threads ; k++ ) {
        r->threads[k].id = k;
        r->threads[k].r = r;
        if ( pthread_create( &r->threads[k].thread , NULL , &runner_main , &r->threads[k] ) != 0 )
            error( "Failed to create runner thread." );
        #if defined(HAVE_SETAFFINITY)
            /* Set the cpu mask to zero | r->id. */
            CPU_ZERO( &cpuset );
            CPU_SET( r->threads[k].id , &cpuset );

            /* Apply this mask to the runner's pthread. */
            if ( pthread_setaffinity_np( r->threads[k].thread , sizeof(cpu_set_t) , &cpuset ) != 0 )
                error( "Failed to set thread affinity." );
        #endif
        }
        
    /* Wait for the runner threads to be in place. */
    while ( r->barrier_count != r->nr_threads )
        if ( pthread_cond_wait( &r->barrier_cond , &r->barrier_mutex ) != 0 )
            error( "Error while waiting for runner threads to get in place." );
    
    }
    
    
    
