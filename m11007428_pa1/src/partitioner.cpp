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

void Partitioner::partition()
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
	// establish _bList
	for(int i = -_maxPinNum; i <= _maxPinNum; i++){
		pair<map<int, Node*>::iterator,bool> bList_ret;
		bList_ret = _bList[0].insert(pair<int, Node*>(i, new Node(-1)));
		if(bList_ret.second == false)
			cout << "bad inser in bList[0]" << endl;
		bList_ret = _bList[1].insert(pair<int, Node*>(i, new Node(-1)));
		if(bList_ret.second == false)
			cout << "bad inser in bList[1]" << endl;
	}
	initialize();
	//TODO:get the _maxGainCell

}

void Partitioner::initialize(){
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
/*	
	for(int i = 0; i < _cellNum; i++){
		cout << "_cellName: " << _cellArray[i]->getName() << endl;
		cout << "_cellGain: " << _cellArray[i]->getGain() << endl;
		cout << "--------------------------------" << endl;
	}
	*/

	int maxGain = INT_MIN;
	for(int i = 0; i < _cellNum; i++){
		int cellGain = 0;
		if(!_cellArray[i]->getPart()){
			cellGain = _cellArray[i]->getGain();
			Node* tempNode = new Node(_cellName2Id[_cellArray[i]->getName()]);
			tempNode->setNext(_bList[0][cellGain]);
			_bList[0][cellGain]->setPrev(tempNode);
			_bList[0][cellGain] = tempNode;
		}
		else{
			cellGain = _cellArray[i]->getGain();
			Node* tempNode = new Node(_cellName2Id[_cellArray[i]->getName()]);
			tempNode->setNext(_bList[1][cellGain]);
			_bList[1][cellGain]->setPrev(tempNode);
			_bList[1][cellGain] = tempNode;
		}
		maxGain = maxGain > cellGain ? maxGain : cellGain;
	}
	cout << "maxGain: " << maxGain << endl;
	if(_bList[0][maxGain]->getId() != -1)
		_maxGainCell = _bList[0][maxGain];
	else
		_maxGainCell = _bList[1][maxGain];
	// printBList();
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

void Partitioner::printBList(){
	for(int i = -_maxPinNum; i <= _maxPinNum; i++){
		cout << "_bList[0][" << i << "]: ";
		Node* tempNode = _bList[0][i];
		while(tempNode != NULL){
			cout << tempNode->getId() << " ";
			tempNode = tempNode->getNext();
		}
		cout << endl;
		
		cout << "_bList[1][" << i << "]: ";
		tempNode = _bList[1][i];
		while(tempNode != NULL){
			cout << tempNode->getId() << " ";
			tempNode = tempNode->getNext();
		}
		cout << endl;
		cout << "-----------------------" << endl;
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
