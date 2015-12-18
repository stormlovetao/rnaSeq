#!/bin/bash
function f1()
{
	echo $0 $1
}
function f2()
{
	echo "f2 #parameters:" $#
	echo "f2 parameters:" $*
	f1 $1 $2
}

f1
