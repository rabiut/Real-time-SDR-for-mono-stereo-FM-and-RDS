## Real-time SDR for mono/stereo FM and RDS

This is the final submission for group 12. We have completed both mono and stereo with all modes functional. A large majority of the RDS has also been completed as it can currently synchronize and detect the blocks. Only the PI code as well as the Program type information can be consistently extracted correctly but we did not have enough time to fully implement program service information and will currently output a mess of some of the decoded characters in the D block.

To run the project code you must first be in the /src subdirectory and use the following command after compilation:

cat "rawIQsamplesfile" | ./project "mode number" "mono stereo selection" | aplay -c "number of channels" -f S16_LE -r "outputFs"

for example:
cat iq_samples.raw | ./project 0 1 | aplay -c 2 -f S16_LE -r 48000
