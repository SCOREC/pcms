#ifndef QUEUE_VISITED_HPP
#define QUEUE_VISITED_HPP

#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_reduce.hpp>
#define MAX_SIZE_QUEUE 500
#define MAX_SIZE_TRACK 800
using namespace std;
using namespace Omega_h;

class queue {
 private:
  LO queue_array[MAX_SIZE_QUEUE];
  int first = 0, last = -1, count = 0;

 public:
  OMEGA_H_INLINE
  queue() {}

  OMEGA_H_INLINE
  ~queue() {}

  OMEGA_H_INLINE
  void push_back(const int& item);

  OMEGA_H_INLINE
  void pop_front();

  OMEGA_H_INLINE
  int front();

  OMEGA_H_INLINE
  bool isEmpty() const;

  OMEGA_H_INLINE
  bool isFull() const;
};

class track {
 private:
  LO tracking_array[MAX_SIZE_TRACK];
  int first = 0, last = -1, count = 0;

 public:
  OMEGA_H_INLINE
  track() {}

  OMEGA_H_INLINE
  ~track() {}

  OMEGA_H_INLINE
  void push_back(const int& item);

  OMEGA_H_INLINE
  int size();

  OMEGA_H_INLINE
  bool notVisited(const int& item);
};

OMEGA_H_INLINE
void queue::push_back(const int& item) {
  if (count == MAX_SIZE_QUEUE) {
    printf("queue is full %d\n", count);
    return;
  }
  last = (last + 1) % MAX_SIZE_QUEUE;
  queue_array[last] = item;
  count++;
}

OMEGA_H_INLINE
void queue::pop_front() {
  if (count == 0) {
    printf("queue is empty\n");
    return;
  }
  first = (first + 1) % MAX_SIZE_QUEUE;
  count--;
}

OMEGA_H_INLINE
int queue::front() { return queue_array[first]; }

OMEGA_H_INLINE
bool queue::isEmpty() const { return count == 0; }

OMEGA_H_INLINE
bool queue::isFull() const { return count == MAX_SIZE_QUEUE; }

OMEGA_H_INLINE
void track::push_back(const int& item) {
  if (count == MAX_SIZE_TRACK) {
    printf("track is full %d\n", count);
    return;
  }
  last = (last + 1) % MAX_SIZE_TRACK;
  tracking_array[last] = item;
  count++;
}

OMEGA_H_INLINE
bool track::notVisited(const int& item) {
  int id;
  for (int i = 0; i < count; ++i) {
    id = (first + i) % MAX_SIZE_TRACK;
    if (tracking_array[id] == item) {
      return false;
    }
  }
  return true;
}

OMEGA_H_INLINE
int track::size() { return count; }

#endif
