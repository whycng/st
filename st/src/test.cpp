
 







// 遍历文件夹，及txt
/*
    std::vector<std::string> my_file;
    std::string need_extension = ".txt";


    // 测试，遍历文件夹下所有文件夹
    std::vector<std::string> folders;
    std::string path = R"(E:\Proj\vsProj\st_FileSave\数据-原始备份 - 副本)";

    get_sub_folders(path, folders);

    // 打印结果
    for (std::string folder : folders) {
        std::cout << folder << std::endl;
        my_file.clear();
        get_need_file(folder, my_file, need_extension);
        for (int i = 0; i < my_file.size(); i++)
        {
            std::cout << "File " << i + 1 << " is:" << std::endl;
            std::cout << my_file[i] << std::endl;
        }
        if (my_file.size() == 0)
        {
            std::cout << "No file can be found!" << std::endl;
        }
        else
        {
            std::cout << std::endl << "Find " << my_file.size() << " file(s)." << std::endl;
        }
    }

*/





//测试，遍历某文件夹下所有txt文件 
/*
std::string file_path = R"(E:\Proj\vsProj\st_FileSave\数据-原始备份 - 副本\1-富源3-H3)";
std::vector<std::string> my_file;
std::string need_extension = ".txt";
   get_need_file(file_path, my_file, need_extension);
    for (int i = 0; i < my_file.size(); i++)
    {
        std::cout << "File " << i + 1 << " is:" << std::endl;
        std::cout << my_file[i] << std::endl;
    }
    if (my_file.size() == 0)
    {
        std::cout << "No file can be found!" << std::endl;
    }
    else
    {
        std::cout << std::endl << "Find " << my_file.size() << " file(s)." << std::endl;
    } */