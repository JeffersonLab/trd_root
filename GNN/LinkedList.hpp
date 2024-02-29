#include <stdio.h>


/**
 * Linked List implementation
 */
class LinkedList{
public:
  class Node{
  public:
    Node* next;
    Node* prev;
    int id;

    Node(int _id){
      this->id = _id;
      next = nullptr;
      prev = nullptr;
    }
  };
  LinkedList(){
    _size = 0;
    head = new Node(-1);
    tail = head;
    head->next = tail;
    head->prev = tail;
    tail->next = head;
    tail->prev = head;
    resetHead();
  }
  void next(){
    current = current->next;
  }
  void prev(){
    current = current->prev;
  }
  void resetHead(){
    current = head;
  }
  void resetTail(){
    current = tail;
  }
  Node getCurrent(){
    return *current;
  }
  Node getHead(){
    return *head;
  }
  Node getTail(){
    return *tail;
  }
  void insert_list_after(LinkedList* ll){
    if(_size == 0){
      head = ll->head;
      tail = ll->tail;
      _size = ll->_size;
      resetHead();
    } else {
      ll->tail->next = current->next;
      ll->head->prev = current;
      current->next->prev = ll->tail;
      current->next = ll->head;
      tail = head->prev;
      _size += ll->_size;
    }
    current = ll->tail;
  }
  void insert_list_before(LinkedList* ll){
    if(_size == 0){
      head = ll->head;
      tail = ll->tail;
      _size = ll->_size;
      resetHead();
    } else {
      ll->head->prev = current->prev;
      ll->tail->next = current;
      current->prev->next = ll->head;
      current->prev = ll->tail;
      head = tail->next;
      _size += ll->_size;
    }
    current = ll->tail;
  }
  void insert_after(int id){
    Node* n = new Node(id);
    if(_size == 0){
      delete head;
      head = n;
      tail = n;
      head->next = tail;
      head->prev = tail;
      tail->next = head;
      tail->prev = head;
    } else {
      n->next = current->next;
      n->prev = current;
      current->next->prev = n;
      current->next = n;
      tail = head->prev;
    }
    current = n;
    _size++;
  }
  void insert_before(int id){
    Node* n = new Node(id);
    if(_size == 0){
      delete head;
      head = n;
      tail = n;
      head->next = tail;
      head->prev = tail;
      tail->next = head;
      tail->prev = head;
    } else {
      n->prev = current->prev;
      n->next = current;
      current->prev->next = n;
      current->prev = n;
      head = tail->next;
    }
    current = n;
    _size++;
  }
  void remove(){
    if(_size == 0){
      throw 1;
    } else if(_size == 1){
      head = new Node(-1);
      tail = head;
      head->next = tail;
      head->prev = tail;
      tail->next = head;
      tail->prev = head;
      resetHead();
    } else {
      Node* next = current->next;
      current->prev->next = current->next;
      current->next->prev = current->prev;
      if(current == head){
        head = current->next;
      } else if(current == tail){
        tail = current->prev;
      }
      delete current;
      current = next;
    }
    _size--;
  }
  int size(){
    return _size;
  }
private:
  Node* head;
  Node* tail;
  Node* current;
  int _size;
};

