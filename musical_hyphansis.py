#!/usr/bin/env python3

import music21
import sys
import math

if len(sys.argv) != 3:
	print('Usage: ./musical_hyphansis source.xml output.txt')
	exit(1)

# First we need to read the musical score 
score = music21.converter.parse(sys.argv[1])
nvoice = len(score.parts)
nmeasure = len(score.parts[0])
print('There are',nvoice,'voices and',nmeasure-1,'measures in this score.')

fhandle = open(sys.argv[2],'w')
for i in range(0,nmeasure):
	for j in range(0,nvoice):
		# Grab the current measure
		cmeasure = score.parts[j].measure(i)
		# If it's empty we're done
		if cmeasure is None:
			continue
		# The algorithm written here assumes that every quarter length of 
		# a note has the form n*2^(-m) where n is a positive integer and 
		# m a non-negative integer, so that we simply need to find the right 
		# power of two to multiply through by in order to ensure that all 
		# of the note durations are integral. 
		# Some example quarter lengths:
		# Whole note         = 4.0    (4*2^0)
		# Half note          = 2.0    (2*2^0)
		# Quarter note       = 1.0    (1*2^0)
		# Eighth note        = 0.5    (1*2^(-1))
		# 16th note          = 0.25   (1*2^(-2))
		# 32nd note          = 0.125  (1*2^(-3))
		# Dotted eighth node = 0.75   (3*2^(-2))
		# Dotted 16th note   = 0.375  (3*2^(-3))
		n = len(cmeasure.notes)
		lmult = 1.0
		for inote in cmeasure.notes:
			length = inote.duration.quarterLength
			f,x = math.modf(length)
			if f > 0.0:
				while True:
					g,x = math.modf(lmult*length)
					if g > 0.0:
						lmult = 2.0*lmult
						# Sanity check...
						if lmult > 32.0:
							print('Error in note duration calculation with',length)
							exit(1)
					else:
						break
		# Now we need to grab the notes for this measure
		for inote in cmeasure.notes:
			# Get this note's frequency
			f = inote.pitch.frequency
			# We operate using an equal-tempered concert piano tuned so that A4 = 440 Hz 
			piano_key = int(round(12.0*math.log(f/440.0,2) + 49.0))
			# Now get the note's duration as an integer by multiplying the quarter length 
			# by the previously computed norm
			length = int(round(lmult*inote.duration.quarterLength))
			# Finally write it out to the file - we divide the voice count by two to handle 
			# the grand staff used for a piano
			for l in range(0,length):
				fhandle.write(str(i) + '/' + str(j/2) + '/' + str(piano_key) + '\n')
fhandle.close()
