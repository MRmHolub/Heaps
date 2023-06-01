
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#define NOT_DEFINED 9999

typedef struct Binom_node{
	int key;	
	struct Binom_node* child; //the leftmost one
	struct Binom_node* sibling;
	struct Binom_node* parent;
	int degree;
} Binomial_Node;

typedef struct {		
	Binomial_Node *root;	
} Binomial_Heap;
// Must be binomial heap not the node, because you want to keep the pointer same, while the pointer to the first node changes by the time


typedef struct {	
	int* values;
	int size;
	int max_size;
} Binary_Heap;

typedef struct Fibbonacci {
	int key;
	struct Fibbonacci* left; //both directional cyclic list
	struct Fibbonacci* right;
	struct Fibbonacci* child; //random one of the children
	struct Fibbonacci* parent;
	int degree; //number of children
	int mark;
} Fibb_Node;

typedef struct {
	Fibb_Node* min;
	int size; //how many nodes are there in heap (not only the main linked list)
}Fibb_Heap;

Binary_Heap* make_binary_Heap(int size) {
	Binary_Heap* h = (Binary_Heap*)malloc(sizeof(Binary_Heap));
	if (h) {			
		h->values = (int*)malloc(size * sizeof(int));
		if (h->values) {
			h->max_size = size;
			h->size = 0;			
			return h;
		}
		free(h);
	}
	printf("Allocation problem !\n\n");
	return NULL;
}


Binomial_Heap* make_binomial_Heap() {
	Binomial_Heap* h = (Binomial_Heap*)malloc(sizeof(Binomial_Heap));
	if (h) {
		h->root = NULL;
		return h;
	}
	printf("Allocation problem !\n\n");
	return NULL;
}
Binomial_Node* make_binomial_Node() {
	Binomial_Node* x = (Binomial_Node*)malloc(sizeof(Binomial_Node));
	if (x) {
		x->sibling = NULL;
		x->parent = NULL;
		x->child = NULL;
		x->degree = 0;
		return x;
	}
	printf("Allocation problem !\n\n");
	return NULL;
}


int parent(int i) {
	return (int)ceil(i / 2.0) - 1;
}

int left(int i) {
	return 2*i + 1;
}

int right(int i) {
	return 2*i + 2;
}

void swap(int* a, int* b) {
	int tmp = *a;
	*a = *b;
	*b = tmp;
}

void min_heapify(Binary_Heap* q, int i){
	/*Vlozil jsem do korene ten nejvetsi prvek a ted ho musim probublat zas dolu, ptam se v kazde vrstve ktery z potomku je mensi a stim swapnu, musim zachovat min haldu*/
	int l = left(i);
	int r = right(i);
	int index_of_smallest;


	if ((l < q->size) && (q->values[l] < q->values[i])) { //if l is not out of the range of my values and is smaller than the root
		index_of_smallest = l;
	}
	else {
		index_of_smallest = i;
	}

	if ((r < q->size) && (q->values[r] < q->values[index_of_smallest])) {
		index_of_smallest = r;
	}

	if (index_of_smallest != i) {
		swap(&q->values[index_of_smallest], &q->values[i]);
		min_heapify(q, index_of_smallest);
	}
}

int binary_remove_min(Binary_Heap* q){
	if (q->size < 1) {
		printf("The heap is already empty, cannot remove min\n");
		return NOT_DEFINED;
	}
	int min = q->values[0];	
	q->values[0] = q->values[q->size - 1]; //No need to delete the last value since I only work with those that are accessible by heap size
	q->size -= 1;
	min_heapify(q, 0);
	
	return min;	
}

void heap_decrease_key(Binary_Heap* q, int key) {
	/*Prvek vlozim nakonec a probublavam nahoru dokud zase neplati ze je halda minimalizacni (kazdy prvek ma pod sebou jen vetsi) */
	int i = q->size-1;
	q->values[i] = key; 
	int tmp;

	while ((i > 0) && (q->values[parent(i)] > q->values[i])){
		tmp = q->values[i];
		q->values[i] = q->values[parent(i)];
		q->values[parent(i)] = tmp;
		i = parent(i);
	}
}

void binary_insert(Binary_Heap* Heap, int key) {
	/*
	We will insert the value at the end and then bubble it upwards until the heap is not correct
	*/
	if (Heap->size >= Heap->max_size) {
		//printf("The Heap is full of: %d, adjusting\n", Heap->max_size);
		int* new_p = realloc(Heap->values, 2 * Heap->max_size * sizeof(int));		
		if (new_p) Heap->values = new_p;
		Heap->max_size = 2 * Heap->max_size;
	} 

	Heap->size = Heap->size + 1;		
	heap_decrease_key(Heap, key);					
}

Binary_Heap* binary_union(Binary_Heap* first, Binary_Heap* second) 
{
	Binary_Heap* h = make_binary_Heap(first->size + second->size);

	for (int i = 0;i < first->size; i++) {
		binary_insert(h, first->values[i]);
	}

	for (int j = 0; j < second->size; j++) {
		binary_insert(h, second->values[j]);
	}	
	return h;
}

void print_binary(Binary_Heap* q) {
	for (int i = 0; i < q->size; i++) {
		printf("%d ", q->values[i]);
	}
	printf(" /n");
}

void free_binary(Binary_Heap* h) {
	if (h) {
		if (h->values) {
			free(h->values);
		}
		free(h);
	}
}


//The binomial source for //extract min, union, insert operations


void merge_binomial_heaps(Binomial_Heap* h1, Binomial_Heap* h2) {	
	//Returns incorrect binomial heap merged in to the first one (incorrect means, there can be 2 trees/siblings of the same order)	
	//Both must contain at least one value

	Binomial_Node* n1, *n2, *to_return;	
	if (h1->root->degree <= h2->root->degree) {	//If the second heap starts with smaller order element we need to swap them (so we dont have to make new Node)
		n1 = h1->root;
		n2 = h2->root;
		to_return = n1;
	} else {
		n1 = h2->root;
		n2 = h1->root;
		to_return = n1;
	}

	Binomial_Node* rest = NULL;
	Binomial_Node* prev = NULL;			

	while (n2 != NULL) {					
		if (!n1) {
			prev->sibling = n2; //I only get there when prev is no longer NULL
			n2 = NULL;
		}
		else if (n1->degree <= n2->degree) {				
			prev = n1;
			n1 = n1->sibling;
		}
		else {
			rest = n2->sibling; //The rest of what I'll need to merge to the node1
			n2->sibling = n1; //I put the rest of the original node to the new one, before inserting the new node to the original
			prev->sibling = n2; //insertion to original
				
			//preparing for the next round
			prev = n2; 
			n2 = rest;
		}
	}
	h1->root = to_return;	
	free(h2);
}

void binomial_link(Binomial_Node* y, Binomial_Node* z) {
	//Connects y node to z as a child node
	y->parent = z;
	y->sibling = z->child;
	z->child = y;
	z->degree += 1;
}

void binomial_union(Binomial_Heap* first, Binomial_Heap* second) {
	//The resulting union will be in the left Heap, the right heap will be freed from memmory	
	if (second && second->root) {
		if (!first->root) {
			first->root = second->root;
		}
		else {
			merge_binomial_heaps(first, second);

			Binomial_Node* prev = NULL;
			Binomial_Node* curr = first->root;
			Binomial_Node* next = curr->sibling;
			Binomial_Node* next_sibling;

			while (next) {
				next_sibling = next->sibling;
				/*I just continue on controling, the order is fine or because there are 3 trees same height (curr, next, nextsib) so I let this one alone, and by merging the other trees I will get a tree one level higher*/
				if ((curr->degree != next->degree) || (next_sibling &&
					(next_sibling->degree == curr->degree))) {
					prev = curr;
					curr = next;
				}
				else {
					if (curr->key <= next->key) { //lets compare which tree will remain and which will be subconnected
						curr->sibling = next->sibling;
						binomial_link(next, curr);
						//We are intentionally not changing prev and stuff
					}
					else {//Subconnecting current node under the next one								
						if (prev == NULL) {
							first->root = next; //This would change the pointer to the heap
						}
						else {
							prev->sibling = next;
						}
						//The actuall connecting
						binomial_link(curr, next);
						//prev = curr; //I am putting to prev the one I put under
						curr = next;
					}
				}
				next = curr->sibling;
			}
		}
	}
}

void binomial_insert(Binomial_Heap* h, int key) {
	//Makes a new bin-tree from key and merges it together
	if (!h->root) {
		h->root = make_binomial_Node();	
	}
	if (h->root->degree == 0) {
		h->root->degree = 1;
		h->root->key = key;
	}
	else {
		Binomial_Heap* neww = make_binomial_Heap();
		neww->root = make_binomial_Node();
		if (neww && neww->root) {
			neww->root->key = key;
			neww->root->degree += 1;
			binomial_union(h, neww);			
		}
	}
}

Binomial_Heap* reverse_binomial_heap(Binomial_Node* root) {	
	Binomial_Node* child = root->child;
	Binomial_Node* curr_child = child;
	Binomial_Node* prev = NULL;
	Binomial_Node* next = NULL;	
	Binomial_Heap* h = make_binomial_Heap();
	//I just change the direction of pointers and for each child in depth 1 delete parent
	if (!h) {
		printf("Not allocated, min not problem");
		return NULL;
	}
	free(root);
	if (curr_child) { //If it is only one node or the memory inicialization goes wrong return NULL
		
		/*Just reversing the pointers */
		while (curr_child)
		{
			curr_child->parent = NULL;//disconnects from the node we are removing

			next = curr_child->sibling; 
			curr_child->sibling = prev; //changes the direction of pointers
			
			prev = curr_child;			
			curr_child = next;														
		}
		h->root = prev;
		return h;
	}
	
	return NULL;
}

int binomial_remove_min(Binomial_Heap* h) {
	/* 
	I will find min value in the root or its siblings (because it is min-heap --> children are always bigger then parents)
	then disconnect it from tree, now the tree without root is actually just reversed binomial heap so you disconnect it, reverse the heap and then union with the original
	*/
	if (!h->root) return NOT_DEFINED;
	Binomial_Node* min = h->root;
	Binomial_Node* before_min_node = NULL;
	Binomial_Node* tmp = h->root;
	
	while (tmp->sibling) {		
		if (tmp->sibling->key < min->key) {
			min = tmp->sibling;
			before_min_node = tmp;
		}	
		tmp = tmp->sibling;
	}

	if (before_min_node) {
		before_min_node->sibling = min->sibling;		//Disconnect from heap if there is a node to disconnect it from, else remove it
	} else {
		h->root = min->sibling;
	}
	int key = min->key;
	binomial_union(h, reverse_binomial_heap(min)); //I am making union on h with removed min and the subheap created be deleting the root (min)
	
	return key;
}

void free_binomial_node(Binomial_Node* h) {
	if (h->child) {
		free_binomial_node(h->child);
	}
	if (h->sibling) {
		free_binomial_node(h->sibling);
	}
	free(h);
}

void free_binomial(Binomial_Heap* h) {
	if (h) {
		if (h->root) {
			free_binomial_node(h->root);
		}
		free(h);
	}
}

void free_fibb_node(Fibb_Node* h) {
	if (h->child) {
		free_fibb_node(h->child);
	}
	if (h->right != h) {	
		if (h->right == h->left) {
			h->right->right = h->right;
			h->right->left = h->right;			
		}
		else {
			h->left->right = h->right;
			h->right->left = h->left;
		}
		free_fibb_node(h->right);
	}
	free(h);
}

void free_fibb(Fibb_Heap* h) {
	if (h) {
		if (h->min) {
			free_fibb_node(h->min);
		}
		free(h);
	}
}

//FIBBONACCI OPERATIONS
//you can have multiple same level trees
//Repairing previous point is delayed until neccessary

Fibb_Heap* make_fibb_Heap() {
	Fibb_Heap* x = (Fibb_Heap*)malloc(sizeof(Fibb_Heap));
	if (x) {
		x->size = 0;
		x->min = NULL;
		return x;
	}
	return NULL;
}

Fibb_Node* make_fibb_node(int key){
	Fibb_Node* node = (Fibb_Node*)malloc(sizeof(Fibb_Node));
	if (node) {
		node->child = NULL;
		node->degree = 0;
		node->mark = 0;
		node->parent = NULL;
		node->left = node;
		node->right = node;
		node->key = key;
		return node;
	}
	return NULL;
}


void fibb_insert(Fibb_Heap* h, int key) {
	//Connect to the both directional list (for example to the left). If the key is the lowest, I must connect heap to it.
	Fibb_Node* node = make_fibb_node(key);
	if (node) {
		if (!h->min) h->min = node;
		else {
			Fibb_Node* root = h->min;
			Fibb_Node* comrade = h->min->left;
			if (h->size == 1) { //both are pointing to each other from both directions				
				root->right = node;
				node->left = root;				
			}
			else {				
				node->left = comrade;				
				comrade->right = node;
			}
			root->left = node; //I am connecting the new node from the left side of the min value in list
			node->right = root;

			if (key < root->key) {
				h->min = node;
			}
		}
		h->size += 1;
	}
}

void conn2v2(Fibb_Node * min1, Fibb_Node * min2) {
	Fibb_Node* min2_end = min2->left;
	Fibb_Node* min1_end = min1->right;
	min1->right = min2;
	min2_end->right = min1_end;
	min2->left = min1;
	min1->left->left = min2_end;
}

void conn1v2(Fibb_Node* min1, Fibb_Node* min2) {
	Fibb_Node* min2_end = min2->left;
	Fibb_Node* min1_end = min1->right;
	min1->right = min2;
	min2_end->right = min1_end;
	min2->left = min1;
	min1->left = min2_end;
}

void conn1v3(Fibb_Node* min1, Fibb_Node* min2) {
	
	Fibb_Node* min2_end = min2->left;
	min1->right = min2;
	min2->left->right = min1;
	min2->left = min1;
	min1->left = min2_end;
}

void conn3v3(Fibb_Node* min1, Fibb_Node* min2) {	
	Fibb_Node* min1_end = min1->right;
	Fibb_Node* min2_end = min2->left;
	min1->right = min2;
	min2->left->right = min1_end;
	min2->left = min1;
	min1_end->left = min2_end;
}

int is_connect(Fibb_Node* min1) {
	Fibb_Node* tmp1 = min1->right;
	
	//printf("\n\n%d %d\n", min1->key, min2->key);
	while (tmp1 != min1) {
		tmp1 = tmp1->right;
		tmp1->mark = 0;
	}
}

void connect(Fibb_Node* min1, Fibb_Node* min2) {
	Fibb_Node* min1_end = min1->right;
	Fibb_Node* min2_end = min2->left;
	Fibb_Node* tmp;

	//does not work for ? Must check for all permutations of 3's, 2's, 1's, 1v2, 1v3 does not work

	if (min1_end == min1 && min2 != min2_end && (min2->right == min2_end)) conn1v2(min1, min2); //1v2
	else if (min1_end != min1 && min2 == min2_end && (min1->left == min1_end)) conn1v2(min2, min1); //2v1
	else if (min1_end == min1 && min2_end != min2->right) conn1v3(min1, min2); //1v3
	else if (min2_end == min2 && min1_end != min1->left) conn1v3(min2, min1); //3v1
	else if (min2_end != min2->right && min1_end == min1->left) conn3v3(min2, min1); //2v3
	else if (min2_end == min2->right && min1_end == min1->left && (min2 != min2_end) && (min1 != min1_end)) conn2v2(min1, min2); //2v2  
	else conn3v3(min1, min2); //it actually works for 3v3, 3v1, 3v2, 1v3, 1V1, 
}

void fibb_union(Fibb_Heap* h1, Fibb_Heap* h2) {
	//Conect h2 into the list of h1; h1 cannot be NULL
	
	Fibb_Node* min1 = h1->min;
	Fibb_Node* min2 = h2->min;
	
	if (!min1) h1->min = min2;
	if (min2) {

		connect(min1, min2);

		if (min1->key > min2->key) { //I can reconnect to the random node since i will get to each of them
			h1->min = min2;
		}
		h1->size += h2->size;
	}
	free(h2);
}

void heap_link(Fibb_Node* y, Fibb_Node* x) {
	//Connect y as a child of x
	//Disconect y from h

	if (y->left == y->right || x->left == x->right) {
		x->left = x;
		x->right = x;
	}
	else {
		y->left->right = y->right; 	//delete y from h, b
		y->right->left = y->left;
	}
	y->parent = x;
	y->left = y;
	y->right = y;

	if (!x->child) {
		x->child = y;		
	}
	else {
		connect(x->child, y);
	}

	if (x->left == y || x->right == y) {
		x->left = x;
		x->right = x;
	}

	x->degree +=1;	
	y->mark = 0;
}

void consolidate(Fibb_Heap* h) {
	//We make sure there is max 1 tree of given order
	
	//CLean up merging
	int max_degree = ceil(log2(h->size));
	Fibb_Node** degrees = calloc(max_degree+2, sizeof(Fibb_Node*)); //Max degree of a node in list + 1
	Fibb_Node* stop = h->min;
	Fibb_Node* curr = h->min;
	Fibb_Node* next = NULL;
	Fibb_Node* same_degree, *tmp;
	Fibb_Node* smallest = h->min;

	int degree;
	do {
		curr->mark = 1;
		next = curr->right; //If only two remain, next is bothways me
		degree = curr->degree;
		while (degrees[degree]) {

			same_degree = degrees[degree];

			if (curr->key > same_degree->key) {//I want to have the small on the top when doing binomial merge, so I must swap them
				tmp = curr;
				curr = same_degree;
				same_degree = tmp;
			}	
			
			heap_link(same_degree, curr);

			degrees[degree] = NULL; //Because by merging two of this degree I get one of the higher degree
			degree += 1;

		}
		degrees[degree] = curr;
		curr = next;
		
		if (curr->key < smallest->key) smallest = curr; // I must find the next min and*
	} while (curr != stop && curr->mark ==0);
	
	is_connect(curr);

	free(degrees);
	h->min = smallest;
}

void point_to_self(Fibb_Node* node) {
	node->left = node;
	node->right = node;
}

int fibb_remove_min(Fibb_Heap* h){
	//In removal the structure is recalculating
	//must take care of calling empty shit in union

	if (h->size) {
		Fibb_Node* z = h->min;

		if (z->degree) { //if the node has children
			Fibb_Heap* children = make_fibb_Heap();
			children->min = z->child;
			children->size = 0;
			Fibb_Node* curr = z->child;

			do { //remove pointers on parent
				curr->parent = NULL;
				curr = curr->right;
			} while (curr != z->child);

			fibb_union(h, children); //I connect children of the min to the list because I will remove the whole tree in next step 		
		}
		int key = z->key;
		if (z->right == z) h->min = NULL; //I removed the last node (there were no children) and no more trees
		else if (z->right == z->left) {
			h->min = z->right;
			point_to_self(z->right);
		}
		else {
			h->min = z->right; //Make heap min some random other value		
			//Disconnect the removed from others in list
			z->left->right = z->right;
			z->right->left = z->left;			
			consolidate(h); 		
		}
		
		h->size -= 1;
		free(z);	
		return key;
	}
	else return NOT_DEFINED;
}

int main()
{

	Binary_Heap* binary = make_binary_Heap(500001);
	Binomial_Heap* binomial = make_binomial_Heap();	
	Fibb_Heap* fibb = make_fibb_Heap();	
	int i, key, res1, res2;

	srand(time(NULL)); //Generates different seq of randoms each time 

	//Generate the array
	int* keys = malloc(sizeof(int) * 500005);
	clock_t start, end;
	double time_binary, time_binomial, time_fibb;
	int n[] = { 500000};
	int HEAP_SIZE;
	//TESTING INSERTS
	
	printf("INSERT:      binary  | binomial | fibbonacci\n");

	for (int j = 0; j < 1; j++) {
		HEAP_SIZE = n[j];

		//for (i = 0; i < HEAP_SIZE + 2; i++) { keys[i] = rand(); }
		for (i = 0; i < HEAP_SIZE + 2; i++) { keys[i] = HEAP_SIZE - i; }

		printf("N = %d  |", HEAP_SIZE);

		//Binary insert
		start = clock();
		for (i = 0; i < HEAP_SIZE; i++) {
			key = keys[i];
			binary_insert(binary, key);
		}
		end = clock();
		time_binary = (double)(end - start) / CLOCKS_PER_SEC;

		//Binomial insert
		start = clock();
		for (i = 0; i < HEAP_SIZE; i++) {
			key = keys[i];
			binomial_insert(binomial, key);
		}
		end = clock();
		time_binomial = (double)(end - start) / CLOCKS_PER_SEC;

		//Fibbonacci insert
		start = clock();
		for (i = 0; i < HEAP_SIZE; i++) {
			key = keys[i];
			fibb_insert(fibb, key);
		}
		end = clock();
		time_fibb = (double)(end - start) / CLOCKS_PER_SEC;

		printf("%.4f s | %.4f s | %.4f s\n", time_binary, time_binomial, time_fibb);
		//TESTING REMOVING
	}



	printf("\n\nEXTRACT MIN: binary | binomial | fibbonacci\n");
	int x, y, z;

	for (int j = 0; j < 1; j++) {
		HEAP_SIZE = n[j];
		printf("N = %d  |", HEAP_SIZE);

		//Binary remove min
		start = clock();
		for (i = 0; i < HEAP_SIZE -1; i++) {
			x = binary_remove_min(binary);
		}
		end = clock();
		time_binary = (double)(end - start) / CLOCKS_PER_SEC;		

		//Binomial remove min
		start = clock();
		for (i = 0; i < HEAP_SIZE - 1; i++) {
			y = binomial_remove_min(binomial);
		}
		end = clock();
		time_binomial = (double)(end - start) / CLOCKS_PER_SEC;

		//Fibbonacci remove min
		start = clock();
		for (i = 0; i < HEAP_SIZE; i++) {
			z = fibb_remove_min(fibb);
		}
		end = clock();
		time_fibb = (double)(end - start) / CLOCKS_PER_SEC;
		printf("%.4f s | %.4f s | %.4f s\n", time_binary, time_binomial, time_fibb);
	}

	free_binary(binary);
	free_binomial(binomial);
	free_fibb(fibb);

	Binary_Heap* binary2 = make_binary_Heap(10001);
	Binomial_Heap* binomial2 = make_binomial_Heap();
	Fibb_Heap* fibb2 = make_fibb_Heap();
	Binary_Heap* binary3 = make_binary_Heap(100001);
	Binomial_Heap* binomial3 = make_binomial_Heap();
	Fibb_Heap* fibb3 = make_fibb_Heap();
	Binary_Heap* neww;

	if (binary2 && binary3 && fibb2 && fibb3 && binomial2 && binomial3) {
		printf("\n\nUNION:    binary | binomial | fibbonacci\n");
		for (int j = 0; j < 1; j++) {
			HEAP_SIZE = n[j];

			//for (i = 0; i < HEAP_SIZE; i++) { keys[i] = rand(); }
			for (i = 0; i < HEAP_SIZE + 2; i++) { keys[i] = HEAP_SIZE - i; }
			printf("N = %d  |", HEAP_SIZE);

			for (i = 0; i < HEAP_SIZE; i++) {
				key = keys[i];
				binary_insert(binary2, key);
				binomial_insert(binomial2, key);
				fibb_insert(fibb2, key);
			}

			//for (i = 0; i < HEAP_SIZE; i++) { keys[i] = rand(); }
			for (i = 0; i < HEAP_SIZE + 2; i++) { keys[i] = HEAP_SIZE - i; }

			for (i = 0; i < HEAP_SIZE; i++) {
				key = keys[i];
				binary_insert(binary3, key);
				binomial_insert(binomial3, key);
				fibb_insert(fibb3, key);
			}

			start = clock();
			neww= binary_union(binary2, binary3);
			end = clock();
			time_binary = (double)(end - start) / CLOCKS_PER_SEC;

			start = clock();
			binomial_union(binomial2, binomial3);
			end = clock();
			time_binomial = (double)(end - start) / CLOCKS_PER_SEC;

			start = clock();
			fibb_union(fibb2, fibb3);		
			end = clock();
			time_fibb = (double)(end - start) / CLOCKS_PER_SEC;
			printf("%.4f s | %.4f s | %.4f s\n", time_binary, time_binomial, time_fibb);
		}
		free(neww);
		free_binary(binary2);
		free_binary(binary3);
		free_binomial(binomial2);
		//free_fibb(fibb2);
	}
	else printf("\nAllocation problem\n");
	return 0;
	
}