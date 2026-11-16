#ifndef PCMS_LOCALIZATION_QUEUE_VISITED_HPP
#define PCMS_LOCALIZATION_QUEUE_VISITED_HPP

#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_reduce.hpp>
#include "pcms/configuration.h"

namespace pcms
{

class Queue
{
private:
  Omega_h::LO queue_array[PCMS_INTERSECTION_QUEUE_SIZE];
  int first = 0, last = -1, count = 0;

public:
  OMEGA_H_INLINE
  Queue() {}

  OMEGA_H_INLINE
  ~Queue() {}

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

class Track
{
private:
  Omega_h::LO tracking_array[PCMS_INTERSECTION_TRACK_SIZE];
  int first = 0, last = -1, count = 0;

public:
  OMEGA_H_INLINE
  Track() {}

  OMEGA_H_INLINE
  ~Track() {}

  OMEGA_H_INLINE
  bool push_back(const int& item);

  OMEGA_H_INLINE
  int size();

  OMEGA_H_INLINE
  bool notVisited(const int& item);
};

OMEGA_H_INLINE
void Queue::push_back(const int& item)
{
  if (count == PCMS_INTERSECTION_QUEUE_SIZE) {
    printf("queue is full %d\n", count);
    return;
  }
  last = (last + 1) % PCMS_INTERSECTION_QUEUE_SIZE;
  queue_array[last] = item;
  count++;
}

OMEGA_H_INLINE
void Queue::pop_front()
{
  if (count == 0) {
    printf("queue is empty\n");
    return;
  }
  first = (first + 1) % PCMS_INTERSECTION_QUEUE_SIZE;
  count--;
}

OMEGA_H_INLINE
int Queue::front()
{
  return queue_array[first];
}

OMEGA_H_INLINE
bool Queue::isEmpty() const
{
  return count == 0;
}

OMEGA_H_INLINE
bool Queue::isFull() const
{
  return count == PCMS_INTERSECTION_QUEUE_SIZE;
}

OMEGA_H_INLINE
bool Track::push_back(const int& item)
{
  if (count == PCMS_INTERSECTION_TRACK_SIZE) {
    // Visited buffer is full. Report failure so callers can stop expanding the
    // search; otherwise new vertices can never be marked visited and the BFS
    // re-queues them forever (infinite loop).
    return false;
  }
  last = (last + 1) % PCMS_INTERSECTION_TRACK_SIZE;
  tracking_array[last] = item;
  count++;
  return true;
}

OMEGA_H_INLINE
bool Track::notVisited(const int& item)
{
  int id;
  for (int i = 0; i < count; ++i) {
    id = (first + i) % PCMS_INTERSECTION_TRACK_SIZE;
    if (tracking_array[id] == item) {
      return false;
    }
  }
  return true;
}

OMEGA_H_INLINE
int Track::size()
{
  return count;
}

} // namespace pcms
#endif // PCMS_LOCALIZATION_QUEUE_VISITED_HPP
