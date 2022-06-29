DT Week 2

Read the files

1. f77 code and comments in good style and easy-read. But it still need to be unified.
	MARCO: That's good. The long term idea would be to convert the entire source code in F90. 
2. f90 code is not in good style.
   1. goto and do notation without unify format. some code have indent but some not
   2. some comments do not format well
   	MARCO: What does it means? 
   3. For readability, some blank line are needed.
   	MARCO: Feel free to add them 
   4. Many people contributed to this project with different code style that's why the code looks difficultly.
   5. I have created a format temple by using HCN/dipole3.f90. But I still cannot push it to the github. I need the access to the git.
   	MARCO: I would suggest you to create a new personal GITHUB repository and add me and Jonathan as contributors.
   6. Use the /=, == ... to replace .ne., .eq ...
   7. I need 3-5 days to complete the f90 code unify work.
   	MARCO: OK you can report this next week 
3. Next I will write a standard format manual for the whole project. Divided into f90 and f77 version. Maybe this is the most simple work, it helps me to read the whole doc and better understand.
4. I find some io performance problems. When I go through the whole project, I will start this work.



f90 unify format

1. Regular indent 6
2. comment notation ! in the front of the line
3. Code block use the --- and === to separate the function, subroutine.
4. Evert function and subroutine have one blank line with before and after content
5. if statement in one line need a blank line after it.
6. Not complete

The code before and after formatting.

![image-20220627015250952](C:\Users\22279\AppData\Roaming\Typora\typora-user-images\image-20220627015250952.png)

![image-20220627030117900](C:\Users\22279\AppData\Roaming\Typora\typora-user-images\image-20220627030117900.png)

![image-20220627031838997](C:\Users\22279\AppData\Roaming\Typora\typora-user-images\image-20220627031838997.png)

![image-20220627031905712](C:\Users\22279\AppData\Roaming\Typora\typora-user-images\image-20220627031905712.png)

	MARCO: I cannot see the figure in the text editor, are you using any particular software?
