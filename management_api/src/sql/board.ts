import { pool } from "../config/pg-pool";
import { zap } from "../helpers/sql";
import { BoardType } from "../sdkgen/api-generated";
import {
	selectBoardSql,
	selectBoardsSharedWithUserSql,
	selectUserOwnedBoardsSql,
} from "./elaborated.query";

async function createBoard(
	userId: string,
	boardId: string,
	name: string,
	type: BoardType,
) {
	await zap.insert("boards", {
		id: boardId,
		ownerId: userId,
		name,
		type,
		createdAt: new Date(),
		updatedAt: new Date(),
		isArchived: false,
	});
}

async function updateBoard(boardId: string, name: string) {
	await zap.update(
		"boards",
		{ name, updatedAt: new Date() },
		{ id: boardId },
	);
}

async function archiveBoard(boardId: string) {
	await zap.update(
		"boards",
		{ isArchived: true, updatedAt: new Date() },
		{ id: boardId },
	);
}

async function getBoard(boardId: string) {
	return (await selectBoardSql.run({ boardId }, pool))[0];
}

async function getUserBoards(userId: string) {
	return selectUserOwnedBoardsSql.run({ userId }, pool);
}

async function getBoardsSharedWithUser(userId: string) {
	return selectBoardsSharedWithUserSql.run({ userId }, pool);
}

async function createBoardSharing(
	id: string,
	boardId: string,
	shareWithUserId: string,
) {
	await zap.insert("board_sharings", {
		id,
		boardId,
		userId: shareWithUserId,
		createdAt: new Date(),
		updatedAt: new Date(),
		isArchived: false,
	});
}

export const boards = {
	createBoard,
	getUserBoards,
	getBoardsSharedWithUser,
	createBoardSharing,
	getBoard,
	updateBoard,
	archiveBoard,
};
// get user boards
// get board
// get user shared boards
// create board
// archive board
// delete board
// edit board
