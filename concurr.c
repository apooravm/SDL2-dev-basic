#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#define NUM_THREADS 4

pthread_mutex_t mutex;
int counter = 0;

// Working with shared data and mutexes
void* increment_counter(void* thread_id) {
    long tid = (long)thread_id;

    // for (int i = 0; i < 100000; i++) {
    //     pthread_mutex_lock(&mutex); // Lock the mutex
    //     counter++; // Critical section
    //     pthread_mutex_unlock(&mutex); // Unlock the mutex
    // }

	pthread_mutex_lock(&mutex);
	counter++;
	pthread_mutex_unlock(&mutex);

    printf("Thread %ld finished.\n", tid);
    pthread_exit(NULL);
}

// Function to be executed by each thread
void* print_message(void* thread_id) {
    long tid = (long)thread_id; // Cast the argument
    printf("Thread %ld: Hello, world!\n", tid);
    pthread_exit(NULL); // Exit the thread
}

void main2() {
	pthread_t threads[NUM_THREADS];
	pthread_mutex_init(&mutex, NULL);

	// Create threads
    for (long t = 0; t < NUM_THREADS; t++) {
        pthread_create(&threads[t], NULL, increment_counter, (void*)t);
    }

    // Wait for threads to finish
    for (int t = 0; t < NUM_THREADS; t++) {
        pthread_join(threads[t], NULL);
    }

	printf("Final counter value: %d\n", counter);

    pthread_mutex_destroy(&mutex); // Clean up the mutex
    pthread_exit(NULL);
}

int main() {
	main2();
	return 0;

	// Store the thread identifiers for each thread created
    pthread_t threads[NUM_THREADS];
	// store the return code from the pthread_create
    int rc;

	// Why long?
	// Since it is being casted to void*, Using long ensures that the size of the data 
	// being cast matches the size of a pointer on most systems.
	// void* is 8bytes while int may still be 4bytes
	// Casting int to void* could truncate it?
    long t;

	// Create 4 threads
    for (t = 0; t < NUM_THREADS; t++) {
        printf("Main: creating thread %ld\n", t);

		// &threads[t] - Pointer to thread_t variable that stores the thread Id
		// NULL - thread attributes?
		// print_message - func to be executed by thread
		// t - Id? Casted to void*
        rc = pthread_create(&threads[t], NULL, print_message, (void*)t);

		// Non zero value returned if pthread_create fails
        if (rc) {
            printf("Error: unable to create thread, %d\n", rc);
            exit(-1);
        }
    }

    // Wait for threads to finish
	// How???
    for (t = 0; t < NUM_THREADS; t++) {
		// Joins the target threads to the main thread
		// joins child thread to parent thread
		// Thus the name
        pthread_join(threads[t], NULL);
    }

    printf("Main: program completed. Exiting.\n");
    pthread_exit(NULL); // Exit the main thread
}

