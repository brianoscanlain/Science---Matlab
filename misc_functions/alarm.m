%Custom Alarm function that can be run on a laptop which runs matlab. 
%This can be extremely useful when we don't have access to an alarm
%and jetlagged.
%
%The function loads  inbuilt audio file so all we need is matlab!
%
%Brian Scanlon, Galway, Feb 2012

function alarm(hour,minute)
%Load the matlab-accompanying audio file:
load chirp;
y1 = y; Fs1 = Fs;
%Configure our audioplayer and assign a handle:
p=audioplayer(y,Fs);

for i=1:864000 %Perform a check every second, for a total of 24 hours:
pause(09)
x=System.DateTime.Now.Hour;
x2=System.DateTime.Now.Minute;
      if x==hour && x2>=minute
         disp('ALARM!!!!!!!!!!!!!!!');
         for iii=1:10
             %audioplayer(y1,Fs1,'sync') % The chirp signal finishes before the
                 play(p, [1 (get(p, 'SampleRate') * 3)])
                 pause(3)
         end

     break

     end

end