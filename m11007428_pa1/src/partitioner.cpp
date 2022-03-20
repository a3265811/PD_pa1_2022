#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <map>
#include "cell.h"
#include "net.h"
#include "partitioner.h"
using namespace std;


void Partitioner::parseInput(fstream& inFile)
{
    string str;
    // Set balance factor
    inFile >> str;
    _bFactor = stod(str);

    // Set up whole circuit
    while (inFile >> str) {
        if (str == "NET") {
            string netName, cellName, tmpCellName = "";
            inFile >> netName;
            int netId = _netNum;
            _netArray.push_back(new Net(netName));
            _netName2Id[netName] = netId;
            while (inFile >> cellName) {
                if (cellName == ";") {
                    tmpCellName = "";
                    break;
                }
                else {
                    // a newly seen cell
                    if (_cellName2Id.count(cellName) == 0) {
                        int cellId = _cellNum;
                        _cellArray.push_back(new Cell(cellName, 0, cellId));
                        _cellName2Id[cellName] = cellId;
                        _cellArray[cellId]->addNet(netId);
                        _cellArray[cellId]->incPinNum();
                        _netArray[netId]->addCell(cellId);
                        ++_cellNum;
                        tmpCellName = cellName;
                    }
                    // an existed cell
                    else {
                        if (cellName != tmpCellName) {
                            assert(_cellName2Id.count(cellName) == 1);
                            int cellId = _cellName2Id[cellName];
                            _cellArray[cellId]->addNet(netId);
                            _cellArray[cellId]->incPinNum();
                            _netArray[netId]->addCell(cellId);
                            tmpCellName = cellName;
                        }
                    }
                }
            }
            ++_netNum;
        }
    }
    return;
}

void Partitioner::partition(fstream& outFile)
{
	// calculate balance bound
	int lower_bound = (1 - _bFactor) / 2 * _cellNum;
	int upper_bound = (1 + _bFactor) / 2 * _cellNum;
	// move half of cells into B part
	// find maximum pin number
	// calculate part size
	for(int i = 0; i < _cellNum; i++){
		if(i >= _cellNum/2){
			_cellArray[i]->move();
			_partSize[1]++;
		}
		else
			_partSize[0]++;
		_maxPinNum = _maxPinNum >= _cellArray[i]->getPinNum() ? _maxPinNum : _cellArray[i]->getPinNum();
	}
	// calculate cut size
	for(int i = 0; i < _netNum; i++){
		vector<int> cellList = _netArray[i]->getCellList();
		bool in_A = false;
		bool in_B = false;
		for(int j = 0; j < cellList.size(); j++){
			if(!_cellArray[cellList[j]]->getPart()){
				in_A = true;
				_netArray[i]->incPartCount(0);
			}
			else{
				in_B = true;
				_netArray[i]->incPartCount(1);
			}
		}
		if(in_A && in_B)
			_cutSize++;
	}
	
	bool firstRound = true;
	while(_maxAccGain != 0 || firstRound){
		firstRound = false;
		_iterNum++;
		// initial partitioner
		_bList[0].clear();
		_bList[1].clear();
		_maxGainCell = NULL;
		_accGain = 0;
		_maxAccGain = 0;
		_moveNum = 0;
		_bestMoveNum = 0;
		_unlockNum[0] = _partSize[0];
		_unlockNum[1] = _partSize[1];
		_moveStack.clear();
		// establish initial _bList
		for(int i = -_maxPinNum; i <= _maxPinNum; i++){
			pair<map<int, Node*>::iterator,bool> bList_ret;
			bList_ret = _bList[0].insert(pair<int, Node*>(i, new Node(-1)));
			if(bList_ret.second == false)
				cout << "bad inser in bList[0]" << endl;
			bList_ret = _bList[1].insert(pair<int, Node*>(i, new Node(-1)));
			if(bList_ret.second == false)
				cout << "bad inser in bList[1]" << endl;
		}
		// reset cell lock and gain
		for(int i = 0; i < _cellNum; i++){
			_cellArray[i]->unlock();
			_cellArray[i]->setGain(0);
		}
		initialize();
		int ini_cutSize = _cutSize;
		for(int i = 0; i < _cellNum; i++){
		/*	cout << "now move: " << i << endl;
			cout << "_accGain: " << _accGain << endl;
			cout << "_maxAccGain: " << _maxAccGain << endl;
			cout << "_moveNum: " << _moveNum << endl;
			cout << "_bestMoveNum: " << _bestMoveNum << endl;
			cout << "_cutSize: " << _cutSize << endl;
			cout << "_maxGainCell:" << _maxGainCell << endl;
			cout << "_maxGainCell ID: " << _maxGainCell->getId() << endl;
			cout << "maxGain: " << _cellArray[_maxGainCell->getId()]->getGain() << endl;
			cout << "---------------------------------\n";*/
			updating(lower_bound, upper_bound);
		}
		cout << "_iterNum: " << _iterNum << endl;
		cout << "_accGain: " << _accGain << endl;
		cout << "_maxAccGain: " << _maxAccGain << endl;
		cout << "_moveNum: " << _moveNum << endl;
		cout << "_bestMoveNum: " << _bestMoveNum << endl;
		cout << "initial _cutSize: " << ini_cutSize << endl;
		recording();
		cout << "new _cutSize: " << _cutSize << endl;
		cout << "---------------------------------\n";
	}
}

void Partitioner::initialize(){
	// cell gain initialize
	for(int i = 0; i < _netNum; i++){
		vector<int> cellList = _netArray[i]->getCellList();
		for(int j = 0; j < cellList.size(); j++){
			int PC_A = _netArray[i]->getPartCount(0);
			int PC_B = _netArray[i]->getPartCount(1);
			if(!_cellArray[cellList[j]]->getPart()){
				if(PC_A == 1)
					_cellArray[cellList[j]]->incGain();
				if(PC_B == 0)
					_cellArray[cellList[j]]->decGain();
			}
			else{
				if(PC_B == 1)
					_cellArray[cellList[j]]->incGain();
				if(PC_A == 0)
					_cellArray[cellList[j]]->decGain();
			}
		}
	}
	// insert all cell into _bList
	// find _maxGainCell
	int maxGain = INT_MIN;
	for(int i = 0; i < _cellNum; i++){
		int cellGain = 0;
		if(!_cellArray[i]->getPart()){
			cellGain = _cellArray[i]->getGain();
			Node* tempNode = new Node(i);
			insertNode(tempNode, 0, cellGain);
			_cellArray[i]->setNode(tempNode);
		}
		else{
			cellGain = _cellArray[i]->getGain();
			Node* tempNode = new Node(i);
			insertNode(tempNode, 1, cellGain);
			_cellArray[i]->setNode(tempNode);
		}
		maxGain = maxGain > cellGain ? maxGain : cellGain;
	}
	if(_bList[0][maxGain]->getId() != -1)
		_maxGainCell = _bList[0][maxGain];
	else
		_maxGainCell = _bList[1][maxGain];
}

void Partitioner::updating(int lower_bound, int upper_bound){
	map<int,int> oldAffectCells;
	Cell* mvCell = _cellArray[_maxGainCell->getId()];
	int from,to;
	vector<int> netList = mvCell->getNetList();
	//cout << "mvCell: " << _maxGainCell->getId() << endl;
	// get from and to side
	if(!mvCell->getPart()){
		from = 0;
		to = 1;
	}
	else{
		from = 1;
		to = 0;
	}
	// cell movement and update some attributes
	_moveStack.push_back(_maxGainCell->getId());
	deleteNode(mvCell->getNode(), mvCell->getPart(), mvCell->getGain());
	delete mvCell->getNode();
	mvCell->setNode(NULL);
	mvCell->move();
	mvCell->lock();
	_cutSize -= mvCell->getGain();
	_partSize[from]--;
	_partSize[to]++;
	_accGain += mvCell->getGain();
	_moveNum++;
	if (_maxAccGain < _accGain){
		_maxAccGain = _accGain;
		_bestMoveNum = _moveNum;
	}
	_unlockNum[from]--;
	
	// updating associated cells
	for(int i = 0; i < netList.size(); i++){
		Net* mvNet = _netArray[netList[i]];
		vector<int> cellList = mvNet->getCellList();
		// before move
		if(mvNet->getPartCount(to) == 0){
			for(int j = 0; j < cellList.size(); j++){
				Cell* affectCell = _cellArray[cellList[j]];
				if(!affectCell->getLock()){
					if(oldAffectCells.count(cellList[j]) == 0)
						oldAffectCells[cellList[j]] = affectCell->getGain();
					affectCell->incGain();
				}
			}
		}
		else if(mvNet->getPartCount(to) == 1){
			for(int j = 0; j < cellList.size(); j++){
				Cell* affectCell = _cellArray[cellList[j]];
				if(!affectCell->getLock() && affectCell->getPart() == to){
					if(oldAffectCells.count(cellList[j]) == 0)
						oldAffectCells[cellList[j]] = affectCell->getGain();
					affectCell->decGain();
				}
			}
		}
		// net change for mvCell
		mvNet->decPartCount(from);
		mvNet->incPartCount(to);
		// after move
		if(mvNet->getPartCount(from) == 0){
			for(int j = 0; j < cellList.size(); j++){
				Cell* affectCell = _cellArray[cellList[j]];
				if(!affectCell->getLock()){
					if(oldAffectCells.count(cellList[j]) == 0)
						oldAffectCells[cellList[j]] = affectCell->getGain();
					affectCell->decGain();
				}
			}
		}
		else if(mvNet->getPartCount(from) == 1){
			for(int j = 0; j < cellList.size(); j++){
				Cell* affectCell = _cellArray[cellList[j]];
				if(!affectCell->getLock() && affectCell->getPart() == from){
					if(oldAffectCells.count(cellList[j]) == 0)
						oldAffectCells[cellList[j]] = affectCell->getGain();
					affectCell->incGain();
				}
			}
		}
	}
	
	// Update _bList
	for(map<int, int>::iterator m_iter = oldAffectCells.begin(); m_iter != oldAffectCells.end(); m_iter++){
		if(m_iter->second != _cellArray[m_iter->first]->getGain()){
			Cell* affectCell = _cellArray[m_iter->first];
			Node* affectNode = affectCell->getNode();
			deleteNode(affectNode, affectCell->getPart(), m_iter->second);
			insertNode(affectNode, affectCell->getPart(), affectCell->getGain());
		}
	}
	
	// find _maxGainCell
	int unbalancedPart = -1;
	if(_partSize[0] < lower_bound)
		unbalancedPart = 0;
	else if(_partSize[0] > upper_bound)
		unbalancedPart = 1;
	else if(_partSize[1] < lower_bound)
		unbalancedPart = 1;
	else if(_partSize[1] > upper_bound)
		unbalancedPart = 0;
	int maxGain = INT_MIN;
	for(int i = 0; i < 2; i++){
		if(unbalancedPart == i)
			continue;
		for(map<int, Node*>::iterator m_iter = _bList[i].begin(); m_iter != _bList[i].end(); m_iter++){
			if(m_iter->second->getId() != -1 && maxGain < m_iter->first){
				_maxGainCell = m_iter->second;
				maxGain = m_iter->first;
			}
		}
	}
}

void Partitioner::recording(){
	// move back bad move cell
	for(int i = _bestMoveNum; i < _moveNum; i++)
		_cellArray[_moveStack[i]]->move();

	// calculate part size
	_partSize[0] = 0;
	_partSize[1] = 0;
	for(int i = 0; i < _cellNum; i++){
		if(!_cellArray[i]->getPart()){
			_partSize[0]++;
		}
		else
			_partSize[1]++;
	}

	// calculate cut size
	_cutSize = 0;
	for(int i = 0; i < _netNum; i++){
		vector<int> cellList = _netArray[i]->getCellList();
		bool in_A = false;
		bool in_B = false;
		_netArray[i]->setPartCount(0,0);
		_netArray[i]->setPartCount(1,0);
		for(int j = 0; j < cellList.size(); j++){
			if(!_cellArray[cellList[j]]->getPart()){
				in_A = true;
				_netArray[i]->incPartCount(0);
			}
			else{
				in_B = true;
				_netArray[i]->incPartCount(1);
			}
		}
		if(in_A && in_B)
			_cutSize++;
	}
}

void Partitioner::insertNode(Node* node, int part, int gain){
	Node* head = _bList[part][gain];
	head->setPrev(node);
	node->setNext(head);
	_bList[part][gain] = node;
	node->setPrev(NULL);
}

void Partitioner::deleteNode(Node* node, int part, int gain){
	if(node->getPrev() != NULL){
		node->getPrev()->setNext(node->getNext());
	}
	else{
		_bList[part][gain] = node->getNext();
	}
	if(node->getNext() != NULL){
		node->getNext()->setPrev(node->getPrev());
	}
	else
		cout << "missing dummy point: _bList[" << part << "][" << gain << "]" << endl;
}

void Partitioner::printSummary() const
{
    cout << endl;
    cout << "==================== Summary ====================" << endl;
    cout << " Cutsize: " << _cutSize << endl;
    cout << " Total cell number: " << _cellNum << endl;
    cout << " Total net number:  " << _netNum << endl;
    cout << " Cell Number of partition A: " << _partSize[0] << endl;
    cout << " Cell Number of partition B: " << _partSize[1] << endl;
    cout << "=================================================" << endl;
    cout << endl;
    return;
}

void Partitioner::reportNet() const
{
    cout << "Number of nets: " << _netNum << endl;
    for (size_t i = 0, end_i = _netArray.size(); i < end_i; ++i) {
        cout << setw(8) << _netArray[i]->getName() << ": ";
        vector<int> cellList = _netArray[i]->getCellList();
        for (size_t j = 0, end_j = cellList.size(); j < end_j; ++j) {
            cout << setw(8) << _cellArray[cellList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::reportCell() const
{
    cout << "Number of cells: " << _cellNum << endl;
    for (size_t i = 0, end_i = _cellArray.size(); i < end_i; ++i) {
        cout << setw(8) << _cellArray[i]->getName() << ": ";
        vector<int> netList = _cellArray[i]->getNetList();
        for (size_t j = 0, end_j = netList.size(); j < end_j; ++j) {
            cout << setw(8) << _netArray[netList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::writeResult(fstream& outFile)
{
    stringstream buff;
    buff << _cutSize;
    outFile << "Cutsize = " << buff.str() << '\n';
    buff.str("");
    buff << _partSize[0];
    outFile << "G1 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 0) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    buff.str("");
    buff << _partSize[1];
    outFile << "G2 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 1) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    return;
}

void Partitioner::printBList(fstream& outFile){
	for(int i = -_maxPinNum; i <= _maxPinNum; i++){
		outFile << "_bList[0][" << i << "]: ";
		Node* tempNode = _bList[0][i];
		while(tempNode != NULL){
			outFile << tempNode->getId() << " ";
			tempNode = tempNode->getNext();
		}
		outFile << endl;
		
		outFile << "_bList[1][" << i << "]: ";
		tempNode = _bList[1][i];
		while(tempNode != NULL){
			outFile << tempNode->getId() << " ";
			tempNode = tempNode->getNext();
		}
		outFile << endl;
		outFile << "-----------------------" << endl;
	}
}

void Partitioner::clear()
{
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        delete _cellArray[i];
    }
    for (size_t i = 0, end = _netArray.size(); i < end; ++i) {
        delete _netArray[i];
    }
    return;
}
